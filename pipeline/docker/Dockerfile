FROM legendzdy/module_snakemake:1.0.0
USER root
COPY runtime/ /
RUN chmod -R a+rX /opt/project/software/gedi
ENV PATH=/opt/project/projects/ISO-SLAM/.envs/iso_slam_v4/bin:/opt/project/software/iso-slam-gedi-env/bin:/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
ENV PYTHONDONTWRITEBYTECODE=1
WORKDIR /work
ENTRYPOINT ["/opt/conda/bin/snakemake"]
CMD ["--help"]
