FROM condaforge/miniforge3:latest
WORKDIR /app
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy
RUN echo "source activate vs_env" > ~/.bashrc
ENV PATH=/opt/conda/envs/vs_env/bin:$PATH 
COPY . .
RUN chmod +x virtual_screening_pipeline.py
CMD [ "python3", "virtual_screening_pipeline.py", "--help" ]