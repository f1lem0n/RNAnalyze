# Przewidywanie struktury II rzędowej mRNA w oparciu o dopasowanie wielosekwencyjne (MSA) aminokwasowe i nukleotydowe oraz dopasowanie strukturalne

> Konserwatywność struktury transkryptów wskazuje na jej istotność we wczesnym fałdowaniu syntetyzowanego białka.

## Procedura
1. Translacja mRNA na sekwencję aminokwasową, identyfikacja białka, wybór odpowiadającej struktury z [RCSB PDB](https://www.rcsb.org/).
2. Wyszukiwanie homologów sekwencyjnych ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) lub strukturalnych ([SCOPe](https://scop.berkeley.edu/)).
3. MSA dla homologów: dopasowanie sekwencji ([clustalw](https://www.genome.jp/tools-bin/clustalw))
   lub struktur ([USalign](https://zhanggroup.org/US-align/)).
4. Konwersja MSA na dopasowanie nukleotydowe.
5. Obliczanie informacji wzajemnej $M_{ij}$.
6. Predykcja struktury transkryptu:
   - na podstawie MSA ([RNAalifold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAalifold.cgi))
   - na podstawie sekwencji transkryptu oraz deklarowanych ograniczeń wynikających z obliczonych
   $M_{ij}$ ([RNAfold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi))

## Zastosowanie

- transkrypt, jego translacja z homologami sekwencyjnymi
- transkrypt, jego translacja z homologami strukturalnymi
- transkrypt, jego translacja, struktura z elementami symetrii np. struktura *2POL*

Występowanie elementów symetrii w strukturze białka sugeruje równoczesne występowanie
pewnych ograniczeń w sekwencji transkryptu. Sekwencja transkryptu powinna zapewnić te same czy
podobne warunki przy konstrukcji symetrycznych domen oraz równocześnie brak wzajemnej
interakcji między fragmentami sekwencjami kodującymi te domeny.

## Wnioski

- Porównanie predykcji struktur mRNA uzyskanych z różnych procedur:
homologii sekwencyjne i strukturalne.
- Analiza dopasowań strukturalnych i sekwencyjnych białka z elementami symetrii domen.
Ocena konsekwencji tych dopasowań na poziomie struktury transkryptu.
- Jaki jest poziom dopasowania sekwencji kodujących tą samą strukturę białek w rodzinie SCOPe?
- Poszukać danych nt. występowania spinek RNA w metodach RNAseq.
Wykorzystać te informacje do weryfikacji predykcji struktury mRNA.
- Wyniki predykcji struktury RNA: najniższa energia i najbardziej prawdopodobna.
Czy są takie same, czy jest duża odległość między nimi?
