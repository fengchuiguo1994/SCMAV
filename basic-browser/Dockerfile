FROM centos:centos7

# usernames and paswords
ENV MYSQL_USER 'mysql'
ENV BASICUSER_MYSQL_PASSWORD ''
ENV BASIC_ADMIN_USERNAME 'root'
ENV BASIC_ADMIN_PASSWORD 'Basic123'
ENV BASIC_ADMIN_EMAIL 'admin@jax.org'

ENV BASIC_DIR=/opt/basic
ENV LANG "en_US.UTF-8"
ENV LD_LIBRARY_PATH $BASIC_DIR/_py/lib/python2.7/site-packages/extsds/

# update and install basic utils
# Install zlib, XZ, pcre
RUN yum -y update && yum install -y epel-release wget gcc && \
    yum group install -y "Development Tools" && \
    yum install -y python-devel python-setuptools python-pip \
    automake \
    pcre \
    pcre-devel \
    xz xz-devel \
    bzip2-devel \
    zlib-devel \
    bison \
    git-core \
    supervisor \
    vim \
    boost boost-devel \
    boost-system boost-filesystem boost-thread \
    && yum clean all

RUN pip install --upgrade pip==18.1
RUN pip install virtualenv==16.7.10
RUN pip install pygments

# copy Basic Browser source code
COPY basic.tar.gz scripts/basic_setup.sh /tmp/
RUN cd /tmp && bash basic_setup.sh

# Install zlib, XZ, pcre
RUN yum -y update &&  yum install -y pcre \
    pcre-devel \
    xz xz-devel \
    bzip2-devel \
    zlib-devel \
    bison \
    git-core \
    supervisor \
    vim \
    boost boost-devel \
    boost-system boost-filesystem boost-thread \
    && yum clean all
    

# Install PCRE
WORKDIR /tmp
RUN git clone --depth=1 --branch v4.0.2 https://github.com/swig/swig.git && \
    cd /tmp/swig && \
    wget https://sourceforge.net/projects/pcre/files/pcre/8.38/pcre-8.38.tar.gz/download -O pcre-8.38.tar.gz && \
    bash Tools/pcre-build.sh \
    && cd /tmp/swig && ./autogen.sh && ./configure && make && make install && cd /tmp && rm -rf swig


# mysql installation
RUN wget -q http://repo.mysql.com/mysql-community-release-el7-5.noarch.rpm && \ 
    rpm -ivh mysql-community-release-el7-5.noarch.rpm && \
    yum -y update && \
    yum -y install mysql-server mysql-devel && \
    yum clean all
RUN pip install mysql-python
ADD confs/bind_0.cnf /etc/mysql/conf.d/bind_0.cnf

# create mysql database
COPY scripts/mysql-run.sh /scripts/mysql-run.sh
RUN chmod 700 /scripts/mysql-run.sh
RUN chown -R mysql:mysql /var/lib/mysql /var/run/mysqld
RUN . /scripts/mysql-run.sh mysqld --datadir=/var/lib/mysql --user=$MYSQL_USER
ADD scripts/run.sh /scripts/run.sh
RUN chmod 755 /scripts/*.sh

# create mysql volume
VOLUME ["/var/lib/mysql"]
VOLUME ["/var/run/mysqld"]
RUN mkdir -p /var/run/mysqld && \
    chown -R mysql:mysql /var/run/mysqld /var/lib/mysql

# Run the setup script
WORKDIR $BASIC_DIR
RUN /usr/bin/python2.7 setup.py $BASIC_DIR/extsds-bin \
        $BASIC_DIR/pytools/ --python=$BASIC_DIR/_py
RUN $BASIC_DIR/_py/bin/python -m pip install pymongo==2.3

# install mongodb
ADD confs/mongodb-org.repo /etc/yum.repos.d/mongodb-org.repo
RUN yum -y install mongodb-org
VOLUME ["/data/db"]
COPY scripts/syncdb.sh /scripts/syncdb.sh
RUN chmod +x /scripts/syncdb.sh && sh /scripts/syncdb.sh

# add supervisor conf
RUN pip install supervisor-stdout
ADD confs/supervisord.conf /etc/supervisord.conf
ADD scripts/run_basic.sh /scripts/run_basic.sh

# Fix missing boost lib for BASIC extsds library
RUN ln -s /usr/lib64/libboost_thread-mt.so.1.53.0 /usr/lib64/libboost_thread-mt.so.5 && \
    ln -s /usr/lib64/libboost_date_time-mt.so.1.53.0 /usr/lib64/libboost_date_time-mt.so.5 && \
    ln -s /usr/lib64/libboost_system-mt.so.1.53.0 /usr/lib64/libboost_system-mt.so.5 && \
    ln -s /usr/lib64/libboost_filesystem-mt.so.1.53.0 /usr/lib64/libboost_filesystem-mt.so.5 && \
    ln -s /usr/lib64/libboost_iostreams-mt.so.1.53.0 /usr/lib64/libboost_iostreams-mt.so.5 && \
    ln -s /usr/lib64/libboost_system.so.1.53.0 /usr/lib64/libboost_system.so.1.40.0

# install gis utils
RUN cd $BASIC_DIR/ds && $BASIC_DIR/_py/bin/python setup.py install

# add genome files
RUN rm -rf /opt/basic/data/genome/* 
COPY files.zip /tmp
RUN cd /tmp && unzip files.zip -d /opt/basic/data/genome/

# cleanup
RUN yum clean all && \
    rm -rf /var/lib/apt/lists/* \
    /tmp/* /var/tmp/*


EXPOSE 8000/tcp

ENTRYPOINT ["/scripts/run.sh"]
