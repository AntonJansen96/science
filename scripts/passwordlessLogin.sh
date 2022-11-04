#!/bin/bash

printf "Please enter [someone]@[location] you want to configure (e.g. anton@htcluster):\n"
read user

printf "\nWhen asked where to save key, press enter.\n"
printf "When asked to overwrite, select no (n).\n"
printf "When asked for your password, enter password.\n\n"

cd ~
ssh-keygen -t rsa
ssh ${user} mkdir -p .ssh
cat .ssh/id_rsa.pub | ssh ${user} 'cat >> .ssh/authorized_keys'
ssh ${user} "chmod 700 .ssh; chmod 640 .ssh/authorized_keys"

# specify directroy 

printf "Passwordless login setup complete.\n"
printf "Please test by performing ssh ${user}...\n\n"
