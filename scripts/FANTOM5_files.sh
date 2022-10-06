#!/bin/bash

parallel -a /data_husihuilke/shared/SMB/miung/Promoter_similarities/Data/FANTOM5_files.txt --jobs 10 wget -P /data_husihuilke/shared/SMB/miung/Promoter_similarities/Data/FANTOM5_files/
#Esto sirve para descargar los links que estan en el archivo de texto, con 10 cores en paralelo, y los guarda en el directorio.
# revisar si agregar -n
