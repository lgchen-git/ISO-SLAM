// Browser QA through Chromium's pipe protocol; no server or package install.
import {spawn} from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';

const report=path.resolve(process.argv[2]);
const root=path.resolve('tmp/browser_qa');
for(const name of ['profile','cache','config','runtime'])fs.mkdirSync(path.join(root,name),{recursive:true});
const browser=spawn(path.join(root,'chromium/chrome'),[
  '--headless','--no-sandbox','--no-zygote','--single-process',
  '--disable-dev-shm-usage','--disable-gpu','--disable-breakpad','--disable-crash-reporter',
  '--disable-background-networking','--disable-component-update','--no-first-run','--no-default-browser-check',
  '--user-data-dir='+path.join(root,'profile'),'--remote-debugging-pipe'],
  {cwd:'/snap/chromium/current/usr/lib/chromium-browser',stdio:['ignore','ignore','pipe','pipe','pipe'],env:{...process.env,LD_LIBRARY_PATH:'/snap/core24/current/usr/lib/x86_64-linux-gnu:/snap/chromium/current/usr/lib/x86_64-linux-gnu',TMPDIR:root,
    XDG_CACHE_HOME:path.join(root,'cache'),XDG_CONFIG_HOME:path.join(root,'config'),XDG_RUNTIME_DIR:path.join(root,'runtime')}});
let pending=new Map(),seq=0,buffer='',stderr='';
const errors=[];
for(const stream of [browser.stdio[3],browser.stdio[4]])stream.on('error',error=>{console.error(error.message,stderr.slice(-3000));process.exit(1)});
browser.on('exit',(code,signal)=>{if(pending.size){console.error('Chromium exited',code,signal,stderr.slice(-3000));process.exit(1)}});
browser.stderr.on('data',d=>{stderr+=d.toString()});
browser.stdio[4].on('data',d=>{buffer+=d.toString();let end;
  while((end=buffer.indexOf('\0'))>=0){const raw=buffer.slice(0,end);buffer=buffer.slice(end+1);if(!raw)continue;
    const message=JSON.parse(raw);if(message.id){const callback=pending.get(message.id);if(callback){pending.delete(message.id);message.error?callback.reject(new Error(JSON.stringify(message.error))):callback.resolve(message.result)}}
    else if(message.method==='Runtime.exceptionThrown')errors.push(message.params.exceptionDetails);
  }});
function command(method,params={},sessionId){return new Promise((resolve,reject)=>{const id=++seq;pending.set(id,{resolve,reject});browser.stdio[3].write(JSON.stringify({id,method,params,...(sessionId?{sessionId}:{})})+'\0')})}
const timeout=setTimeout(()=>{console.error('Browser QA timeout',stderr.slice(-3000));browser.kill('SIGKILL');process.exit(1)},55000);
try{
  const {targetId}=await command('Target.createTarget',{url:'about:blank'});
  const {sessionId}=await command('Target.attachToTarget',{targetId,flatten:true});
  await command('Runtime.enable',{},sessionId);await command('Page.enable',{},sessionId);
  await command('Emulation.setDeviceMetricsOverride',{width:1440,height:1000,deviceScaleFactor:1,mobile:false},sessionId);
  await command('Page.navigate',{url:'file://'+report},sessionId);
  await new Promise(resolve=>setTimeout(resolve,5000));
  async function evaluate(expression){const r=await command('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true},sessionId);if(r.exceptionDetails)throw new Error(JSON.stringify(r.exceptionDetails));return r.result.value}
  const initial=await evaluate(`({ready:document.readyState,plotly:!!window.Plotly,sections:document.querySelectorAll('.report-section').length,expand:document.querySelectorAll('[data-report-expand]').length,rows:document.querySelectorAll('#result-body tr').length})`);
  if(initial.ready!=='complete'||!initial.plotly||!initial.expand)throw new Error('Incomplete report initialization '+JSON.stringify({initial,errors}));
  const modal=await evaluate(`(()=>{document.querySelector('[data-report-expand]').click();const opened=document.getElementById('report-figure-modal').classList.contains('open');document.dispatchEvent(new KeyboardEvent('keydown',{key:'Escape'}));return {opened,closed:!document.getElementById('report-figure-modal').classList.contains('open')}})()`);
  if(!modal.opened||!modal.closed)throw new Error('Figure modal failed');
  const print=await evaluate(`(()=>{const before=document.querySelectorAll('#result-body tr').length;window.dispatchEvent(new Event('beforeprint'));const all=document.querySelectorAll('#result-body tr').length;window.dispatchEvent(new Event('afterprint'));return {before,all,restored:document.querySelectorAll('#result-body tr').length}})()`);
  if(initial.rows&&(!print.all||print.restored!==print.before||print.all<print.before))throw new Error('Print row handling failed '+JSON.stringify(print));
  const generic=await evaluate(`(()=>{const p=Array.from(document.querySelectorAll('[data-pager-for]')).find(p=>p.dataset.pagerFor!=='project-result-table');if(!p)return null;const t=document.getElementById(p.dataset.pagerFor),rows=Array.from(t.querySelectorAll('tbody tr'));const before=rows.filter(r=>!r.hidden).length;p.querySelector('[data-page-next]').click();const page=p.querySelector('[data-page-status]').textContent;window.dispatchEvent(new Event('beforeprint'));const printed=rows.filter(r=>!r.hidden).length;window.dispatchEvent(new Event('afterprint'));return {before,total:rows.length,printed,page,restored:rows.filter(r=>!r.hidden).length}})()`);
  if(generic&&(generic.printed!==generic.total||generic.restored>10))throw new Error('Generic table print/pagination failed '+JSON.stringify(generic));
  if(initial.rows){
    const interaction=await evaluate(`(()=>{const next=document.getElementById('next-page');next.click();const page=document.getElementById('page-status').textContent;const select=document.getElementById('layer-select');select.selectedIndex=1;select.dispatchEvent(new Event('change'));const filter=document.getElementById('result-filter');filter.value='G1';filter.dispatchEvent(new Event('input'));const filtered=document.getElementById('page-status').textContent;next.click();return {page,filtered,preserved:filter.value==='G1'}})()`);
    if(!interaction.preserved)throw new Error('Pagination cleared filter');
  }
  if(errors.length)throw new Error('Browser JavaScript exceptions '+JSON.stringify(errors));
  const screenshot=await command('Page.captureScreenshot',{format:'png'},sessionId);
  fs.writeFileSync(path.join(root,path.basename(report)+'.png'),Buffer.from(screenshot.data,'base64'));
  console.log(JSON.stringify({status:'PASS',report,initial,modal,print,generic,js_exceptions:errors.length}));
}catch(error){console.error(error.stack,stderr.slice(-3000));process.exitCode=1}
finally{clearTimeout(timeout);browser.kill('SIGTERM')}
