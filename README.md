# Preklapanje sekvenci dugih greškovitih očitanja

## Buildanje i pokretanje
Projekt bi se trebao moći jednostavno izgraditi pozicioniranjem u korijenski direktorij te izvođenjem jedne od skripti za izgradnju (_build.sh_ ili _build.bat_).
Primjer izgradnje i pripreme projekta:  
```sh
$ scripts/build.sh
```  
Primjer pokretanja iz korijenskog direktorija projekta:
```sh
$ build/bin/mapper src/resources/lambda_reads.fasta out.paf 
```  

Postoji mogućnost ograničenja broja dretvi trećim argumentom. Na primjer, ako se program želi ograničiti na 3 dretve, poziva se ovako:
```sh
$ build/bin/mapper src/resources/ecoli_reads.fasta out.paf 3
```  

Detaljnije upute za uporabu mogu se dobiti uz pozivni parametar `-h`:
```sh
$ mapper -h
```  
