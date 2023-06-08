#!/bin/bash

sudo R CMD REMOVE robets
sudo R CMD build robets/
sudo R CMD INSTALL --no-lock robets_1.4.tar.gz

