#!/bin/sh
#SBATCH --time=08:00:00
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --output=rstudio-server.job.%j
#SBATCH --job-name=RStudio_Session

#Set the container for the session
#IMAGE=/opt/singularity-images/nextflow/denom-pkdp-latest.img
#IMAGE=/opt/singularity-images/nextflow/statgen-Jul19.img
IMAGE=/opt/singularity-images/nextflow/denom-bioconductor_docker_genomics-latest.img

#Set a temp password for the session
#Could probably use LDAP for authenicatino...somehow
export PASSWORD=$(openssl rand -hex 4)

# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')

# Extract the walltime to a variable
export SLURM_WALLTIME=$(squeue -j ${SLURM_JOB_ID} -h -o %L)

#Log some vitals to stdout
echo "Session details:"
echo --------------------------------------------
echo "Host: http://${HOSTNAME}.coh.org:${PORT}"
echo "Username: ${USER} Password: ${PASSWORD}"
echo "Slurm JobID: ${SLURM_JOB_ID}"
echo "Singularity image: ${IMAGE}"


#email connection instructions
cat <<END | /usr/sbin/sendmail -t 
To: $USER@coh.org 
Subject: Your RStudio Session on Apollo


Your RStudio session has been submitted to the cluster and is running. 

1. To connect to the session,

   Point your web browser to http://${HOSTNAME}.coh.org:${PORT}

2. Log in to RStudio Server using the following credentials:

   user: ${USER}
   password: ${PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exiting the RStudio Session ("power" button in the top right corner of the RStudio window)

2. Issuing the following command on the login node:

      scancel -f ${SLURM_JOB_ID}

#######################################################
Cluster Resource Request:
CPUs: $SLURM_JOB_CPUS_PER_NODE | Memory: $SLURM_MEM_PER_NODE MB | Time: $SLURM_WALLTIME remaining
Singularity image: $(basename ${IMAGE})
#######################################################

END

# User-installed R packages go into their home directory
if [ ! -e ${HOME}/.Renviron ]
then
  printf '\nNOTE: creating ~/.Renviron file\n\n'
  echo 'R_LIBS_USER=~/R/%p-library/%v' >> ${HOME}/.Renviron
fi

#load singularity
module load singularity

# This example bind mounts the /project directory on the host into the Singularity container.
# By default the only host file systems mounted within the container are $HOME, /tmp, /proc, /sys, and /dev.
mkdir -p /tmp/$USER && singularity run --bind /tmp/$USER:/tmp $IMAGE  \
    rserver --www-port ${PORT} --auth-none=0 --auth-pam-helper-path=pam-helper 
    
printf 'rserver exited' 1>&2
