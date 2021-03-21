#maintainer Ariella Gladstein
#organization tskit-dev
#application "Msprime: A reimplementation of Hudson's classical ms simulator for modern data sets."

FROM jupyter/scipy-notebook:latest

# Set the working directory to /app
WORKDIR /app

# Install latest msprime release
RUN pip install --pre msprime
