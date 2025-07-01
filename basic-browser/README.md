# Basic Browser

[Refer to the TheJacksonLaboratory basic-browser](https://github.com/TheJacksonLaboratory/basic-browser)

## Running Basic Browser
Below are the instructions to run Basic Browser docker image.

## Install Docker
```
yum install docker # or
yum install docker-ce

systemctl start docker.service
systemctl enable docker.service
```

## How to build the Basic Browser Docker image
```sh
docker build -t basic-browser .
cd ..
```

## How to run Basic Browser
This creates a mongodb and mysql Docker volumes, which is useful for data persistence
```sh
docker run -d --name basic-browser -p 8000:8000 -v mongodb:/data/db -v mysql:/var/lib/mysql  basic-browser # basic browser start
```
The basic browser should be running at localhost:8000. To upload track and data into the Basic Browser, a new TTY terminal must be opened in the container.
```sh
docker exec -ti basic-browser bash  # enter basic browser (in docker)
```
Create a username and password through the terminal the first time
```sh
# this should create the superuser 'root'
cd /opt/basic/
_py/bin/python manage.py createsuperuser # give username, emall and password
exit # exit basic browser (in docker)

docker stop basic-browser # stop basic browser
docker start basic-browser # start basic browser
docker restart basic-browser # restart basic browser

docker rm -f $(docker ps --filter name=basic-browser -q) # delete the basic browser container
```

## Install Genome Annotation on Basic Browser
(1) Make sure Basic Browser is running at localhost:8000, and you are able to login using your username and password, and a TTY terminal is opened in the container.<br/>
(2) To create a folder to save gene annotation, go to admin page -> Librarys -> Add library -> Create a library name (for example, "annotation") -> save<br/>
(3) Copy the file to be uploaded inside basic docker <br/>
```sh
docker exec -ti basic-browser bash
scp usename@ip:uploadfile . # hg19_known_genes.txt
```
(4) Install the genome assembly <br/>
```sh
_py/bin/python console/table_util.py create hg19 "UCSC Known Genes (hg19)"
ID=1
_py/bin/python console/table_util.py load_genes ${ID} -i hg19_known_genes.txt
_py/bin/python console/track_util.py gen_genes hg19
```

## Upload Tracks to Basic Browser
(1) Copy the file to be uploaded inside basic docker <br/>
```sh
docker exec -ti basic-browser bash
scp usename@ip:uploadfile . 
```
(2) To create a folder to organize tracks, go to admin page -> Librarys -> Add library -> Create a library name (for example, "test_upload") -> save <br/>

#### Upload Coverage Track 
(1) To upload the coverage track, first create a table for coverage, select a folder ("test_upload"), and give the name of the track ("GM12878_RNAPII_coverage")
```
alias TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py" 
alias TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
TABLE create hg19 -l "test_upload" "GM12878_RNAPII_coverage"
```
(2) Then you will get a TABLE ID for this track, and save the ID to the corresponding variable:
```
COV=2
```
(3) Upload coverage track:
```
TRACK gen_cov max ${COV} /Documents/GM12878_RNAPII_insitu.for.BROWSER.bedgraph
```

#### Upload Peak Track 
(1) To upload the peak track, first create a table for peak, select a folder ("test_upload"), and give the name of the track ("GM12878_RNAPII_peak")
```
alias TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py" 
alias TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
TABLE create hg19 -l "test_upload" "GM12878_RNAPII_peak"
```
(2) Then you will get a TABLE ID for this track, and save the ID to the corresponding variable:
```
BED=4
```
(3) Upload peak track:
```
TABLE load ${BED} 1:chrom 2:start 3:end -i /Documents/GM12878_RNAPII_insitu.for.BROWSER.spp.z6.broadPeak
TRACK new ${BED} scls
```

#### Upload Loop Track 
```
alias TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py" 
alias TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
TABLE create hg19 -l "test_upload" "GM12878_RNAPII_loop"
CLU=3
TABLE load ${CLU} 1:chrom 2:start 3:end 4:chrom2 5:start2 6:end2 7:score -i /Documents/GM12878_RNAPII_insitu.e500.clusters.cis.BE3
TRACK new ${CLU} curv
```

#### Configure Loop Track
Go to Tracks in Admin page and select the loop track you just uploaded, under Metadatas "Key" and "Value", set: 
```sh
Key="options",
Value=
{
  "opacity": 0.02,
  "outbound": {
    "show": false
  },
"yaxis": {
"log": 10.0 }
}
```
Go to the next row,
```sh
Key="series",
Value= [{"color": "red"}]
```