FROM alpine:3.19.1

WORKDIR /app

# install dependencies
RUN apk add wget git tar gawk grep python3 make g++
RUN git clone https://github.com/pylelab/USalign.git
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN tar -xvzf ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN rm -rf ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN cd USalign && make

# update env variables
ENV PATH="/app/USalign:/app/ncbi-blast-2.15.0+/bin:${PATH}"

# copy source code
COPY sequence_homology /app/sequence_homology
COPY structural_homology /app/structural_homology
COPY run.py /app/run.py

CMD ["python3", "run.py"]