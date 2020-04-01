#maintainer Ariella Gladstein
#organization tskit-dev
#application "Msprime: A reimplementation of Hudson's classical ms simulator for modern data sets."

FROM ubuntu:18.04

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

RUN apt-get update && apt-get install -y --no-install-recommends \
            python3 \
            python3-dev \
            python3-pip \
            libgsl-dev \
            build-essential \
            python3-wheel \
            git \
            && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade setuptools
RUN pip3 install .
