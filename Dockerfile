FROM debian:bullseye-slim

WORKDIR /app

# install dependencies
RUN apt-get update && apt-get install -y wget git tar gawk grep python3 make g++ clustalw
RUN git clone https://github.com/pylelab/USalign.git
RUN cd USalign && make
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN rm ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN mv ncbi-blast-2.15.0+ ncbi-blast

# update env variables
ENV PATH="/app/USalign:/app/ncbi-blast/bin:${PATH}"

# copy source code
COPY sequence_homology /app/sequence_homology
COPY structural_homology /app/structural_homology
COPY run.py /app/run.py

CMD ["python3", "run.py"]