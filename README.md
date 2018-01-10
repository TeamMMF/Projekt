# Preklapanje sekvenci dugih greškovitih očitanja
## Verzije
Spori algoritmi za jednu dretvu, ali najbrže vrijeme izvođenja s više procesora:
```sh
build/bin/mapper_parallel src/resources/ecoli_reads.fasta out_t.paf
```  

Radi isključivo jedna dretva:
```sh
build/bin/mapper_v3 src/resources/ecoli_reads.fasta out_t.paf
```  
Budući da su sve konstante određene uglavnom nabadanjem, ovo se može koristiti ako prve dvije verzije daju neočekivano loše rezultate (jaccard ispod 0.8) za nešto što nije ecoli ili lambda:
```sh
build/bin/mapper_dynamic src/resources/ecoli_reads.fasta out_t.paf
```  

## Buildanje i pokretanje
Projekt bi se trebao moći jednostavno izgraditi pozicioniranjem u korijenski direktorij te izvođenjem jedne od skripti za izgradnju (_build.sh_ ili _build.bat_).
Primjer izgradnje i pripreme projekta:  
```sh
$ scripts/build.sh
```  
Primjer pokretanja iz korijenskog direktorija projekta (default je verzija _mapper_parallel_):
```sh
$ build/bin/mapper src/resources/lambda_reads.fasta out.paf 
```  
Detaljnije upute za uporabu mogu se dobiti uz pozivni parametar `-h`:
```sh
$ mapper -h
```  
