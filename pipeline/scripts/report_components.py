"""Shared offline figure, pagination and print behavior."""
import base64
import re
from pathlib import Path

def enhance_report(content, expandable_figures=True):
    content = content.replace('<section id=', '<section class="report-section" id=')
    if expandable_figures:
        content = content.replace('class="figure"', 'class="figure-card" data-expandable="true"')
        figure_pattern = r'(<(?:article|div)\b[^>]*class="[^"]*\bfigure-card\b[^"]*"[^>]*>)'
        content = re.sub(figure_pattern, lambda m: m.group(1) + '<button type="button" class="figure-expand" data-report-expand="true">放大图形</button>', content)
    # Existing compact tables must expose a real pager when they contain results.
    index = 0
    def table(match):
        nonlocal index
        attrs, body = match.groups()
        if 'data-page-size' not in attrs:
            return match.group(0)
        index += 1
        table_id = f'report-table-{index}'
        return f'<table{attrs} id="{table_id}">{body}</table><div class="pager" data-pager-for="{table_id}"><button type="button" data-page-prev>上一页</button><span data-page-status></span><button type="button" data-page-next>下一页</button></div>'
    content = re.sub(r'<table([^>]*)>(.*?)</table>', table, content, flags=re.S)
    content = content.replace('<table><thead><tr><th>特征</th>', '<table id="project-result-table" class="paginated-table" data-page-size="10"><thead><tr><th>特征</th>')
    content = content.replace('<div class="pager"><button id="prev-page"', '<div class="pager" data-pager-for="project-result-table"><button data-page-prev id="prev-page"')
    content = content.replace('id="page-status"', 'data-page-status id="page-status"').replace('id="next-page"', 'data-page-next id="next-page"')
    logo = Path(__file__).resolve().parents[1]/'assets/epi-logo.png'
    if 'EPIBIOTEK' not in content and logo.exists():
        uri='data:image/png;base64,'+base64.b64encode(logo.read_bytes()).decode()
        content=content.replace('<footer class="footer">', f'<footer class="footer"><img src="{uri}" alt="EPIBIOTEK 麦观生物" style="display:block;width:240px;max-width:90%;margin:0 auto 12px">')
    expand_css = '.figure-expand{display:block;margin:0 0 8px auto;background:#fff;border:1px solid #bdc8d5;border-radius:4px;padding:4px 10px;cursor:zoom-in}.report-figure-modal{display:none;position:fixed;inset:0;background:#000b;z-index:1000;padding:3vh 3vw}.report-figure-modal.open{display:flex;align-items:center;justify-content:center}.report-modal-body{background:#fff;padding:18px;width:94vw;max-height:94vh;overflow:auto}.report-modal-content svg,.report-modal-content img{max-width:100%;height:auto}' if expandable_figures else ''
    print_expand = '.report-figure-modal,.figure-expand,' if expandable_figures else ''
    css = f'''<style>{expand_css}.pager{{display:flex;gap:10px;align-items:center;margin:12px 0}}@media print{{.sidebar,.pager,.modal,{print_expand}.selectors,.modebar{{display:none!important}}main,.footer{{margin-left:0!important}}details{{display:block}}section,.figure-card{{break-inside:avoid}}tbody tr{{display:table-row!important}}.figure-card{{box-shadow:none}}.plot{{display:none!important}}.static-fallback{{display:block!important}}}}</style>'''
    extra = '''<div id="report-figure-modal" class="report-figure-modal" role="dialog" aria-modal="true" aria-label="放大图形"><div class="report-modal-body"><button type="button" id="report-modal-close">关闭</button><div class="report-modal-content"></div></div></div><script>
(()=>{const modal=document.getElementById('report-figure-modal'),area=modal.querySelector('.report-modal-content');let focus;
function close(){modal.classList.remove('open');if(focus)focus.focus()}
document.getElementById('report-modal-close').onclick=close;modal.onclick=e=>{if(e.target===modal)close()};document.addEventListener('keydown',e=>{if(e.key==='Escape')close()});
document.addEventListener('click',e=>{const button=e.target.closest('[data-report-expand]');if(!button)return;focus=button;const card=button.closest('.figure-card'),clone=card.cloneNode(true);clone.querySelectorAll('button').forEach(x=>x.remove());clone.querySelectorAll('[id]').forEach(x=>x.removeAttribute('id'));area.replaceChildren(clone);modal.classList.add('open');document.getElementById('report-modal-close').focus();
const source=card.querySelector('.js-plotly-plot'),target=clone.querySelector('.js-plotly-plot');if(source&&target&&window.Plotly){target.replaceChildren();Plotly.newPlot(target,source.data,source.layout,{responsive:true,displaylogo:false})}});
document.querySelectorAll('[data-pager-for]').forEach(pager=>{const table=document.getElementById(pager.dataset.pagerFor);if(!table||table.id==='project-result-table')return;const rows=Array.from(table.querySelectorAll('tbody tr'));let page=1;const size=10;function render(){const pages=Math.max(1,Math.ceil(rows.length/size));page=Math.min(page,pages);rows.forEach((r,i)=>r.hidden=i<(page-1)*size||i>=page*size);pager.querySelector('[data-page-status]').textContent=`第 ${page} / ${pages} 页，共 ${rows.length} 条`;pager.querySelector('[data-page-prev]').disabled=page===1;pager.querySelector('[data-page-next]').disabled=page===pages}pager.querySelector('[data-page-prev]').onclick=()=>{page--;render()};pager.querySelector('[data-page-next]').onclick=()=>{page++;render()};render()});
window.addEventListener('beforeprint',()=>{document.querySelectorAll('details').forEach(d=>{d.dataset.wasOpen=d.open;d.open=true});document.querySelectorAll('tbody tr').forEach(r=>{r.dataset.wasHidden=r.hidden;r.hidden=false})});window.addEventListener('afterprint',()=>{document.querySelectorAll('details').forEach(d=>d.open=d.dataset.wasOpen==='true');document.querySelectorAll('tbody tr').forEach(r=>r.hidden=r.dataset.wasHidden==='true')});
})();</script>'''
    if not expandable_figures:
        extra = ''
    return content.replace('</head>', css+'</head>').replace('</body>',extra+'</body>')
