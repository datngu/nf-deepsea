## BUILD DOCKER AUTOMATICLY
docker build -t deepsea:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag deepsea:v0.0.0 ndatth/deepsea:v0.0.0
docker push ndatth/deepsea:v0.0.0
echo DONE


### test docker

docker run -it --rm -v /sigma4:/sigma4 --name deepsea ndatth/deepsea:v0.0.0
docker start deepsea
docker attach deepsea