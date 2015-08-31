#!/bin/bash

##-- okay, this boot script gets called and run on all worker nodes
##derived from install_rmr2.sh from Jamie Olson
## by Nick Bruns, on Monday, June 21, 2015


# The R version in the ubuntu repo is way too old
# So add a CRAN repo
sudo su << EOF1 
echo ' 
deb http://cran.rstudio.com/bin/linux/ubuntu precise/
' > /etc/apt/sources.list.d/r-cran.list 
EOF1

# If running trusty, you also need to add this to get around
# a missing package
# sudo add-apt-repository ppa:marutter/rrutter

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo apt-get update


branch=master
sudo apt-get install -y r-base-core
#sudo apt-get install -y r-cran-rcpp # Don't install this from the repo, use CRAN directly

sudo R --no-save << EOF
install.packages(c('pryr','PresenceAbsence', 'verification', 'gbm', 'scam', 'splines'), repos="http://cran.revolutionanalytics.com", INSTALL_opts=c('--byte-compile') )
EOF


##, instead of the rmr packages, going to first install my own packages



# Setup environmental variables
sudo su < EOF1 
echo ' 
# Should probably figure out where the hdinsight HADOOP_HOME should be
#export HADOOP_HOME=/usr/lib/hadoop
export HADOOP_CMD=/usr/bin/hadoop
export HADOOP_STREAMING=/usr/hdp/current/hadoop-mapreduce-client/hadoop-streaming.jar
' >> /etc/profile.d/rhadoop.sh 
EOF1

