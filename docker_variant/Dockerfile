FROM condaforge/miniforge3:latest
WORKDIR /app
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean -afy
ENV PATH=/opt/conda/envs/vistadock_p_env/bin:$PATH 
RUN wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /usr/local/bin/gnina && chmod +x /usr/local/bin/gnina
COPY . .
RUN chmod +x vistadockp.py
CMD [ "python3", "vistadockp.py", "--help" ]
