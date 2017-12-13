# Preklapanje sekvenci dugih greškovitih očitanja
Glavni program nalazi se u datoteci _mapper.cpp_.
Projekt bi se trebao moći jednostavno izgraditi pozicioniranjem u korijenski direktorij te izvođenjem jedne od skripti za izgradnju (_build.sh_ ili _build.bat_).
Primjer izgradnje i pripreme projekta:  
```sh
$ scripts/build.sh
```  
Primjer pokretanja iz korijenskog direktorija projekta:
```sh
$ build/bin/mapper src/resources/lambda_reads.fasta out.paf 
```  
Detaljnije upute za uporabu mogu se dobiti uz pozivni parametar `-h`:
```sh
$ mapper -h
```  
