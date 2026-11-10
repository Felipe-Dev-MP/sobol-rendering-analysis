\# Análise de amostragem Sobol 2D



\## Teste de 64 pontos



O teste verifica se a Sobol 2D convencional (`i2x2y`) gera a mesma sequência

que a opção `pbrt` do programa do autor. Também verifica a quantidade, a

unicidade e o intervalo das coordenadas de todos os métodos testados.



\### Como repetir



Baixe `paper1031\_supp.zip` no \[suplemento oficial do artigo](https://diglib.eg.org/items/57f2cdeb-69d9-434e-8cf8-37b63e7e69d9).

Extraia o arquivo e deixe `benchmarking.cpp` e `sobolmatrices.h` na mesma pasta.

Abra o terminal \*\*MSYS2 MSYS\*\* nessa pasta. Se precisar instalar o compilador,

execute `pacman -S --needed gcc`.



Compile o programa e gere os pontos:



```bash

g++ -std=gnu++17 -O2 -include ctime benchmarking.cpp -o benchmarking.exe

mkdir -p resultados



for metodo in i2x2y lp mmd pbrt; do

&#x20; ./benchmarking.exe -a "$metodo" -m 6 "resultados/$metodo.txt" || break

done

```



O `-m 6` pede 2⁶ pontos. A primeira linha de cada arquivo informa a quantidade

de pontos. Cada linha seguinte contém um par de coordenadas.



Para conferir a quantidade e as coordenadas, e comparar os dois arquivos Sobol:



```bash

awk '

&#x20; FNR == 1 { if ($1 != 64) erro = 1; next }

&#x20; { total\[FILENAME]++ }

&#x20; NF != 2 || $1 < 0 || $1 >= 1 || $2 < 0 || $2 >= 1 { erro = 1 }

&#x20; !vistos\[FILENAME, $1 + 0, $2 + 0]++ { distintos\[FILENAME]++ }

&#x20; END {

&#x20;   for (i = 1; i < ARGC; i++) {

&#x20;     arquivo = ARGV\[i]

&#x20;     printf "%s: %d pontos, %d distintos\\n", arquivo, total\[arquivo], distintos\[arquivo]

&#x20;     if (total\[arquivo] != 64 || distintos\[arquivo] != 64) erro = 1

&#x20;   }

&#x20;   print erro ? "Falha na validacao." : "Todos os pontos estao no intervalo esperado."

&#x20;   exit erro

&#x20; }

' resultados/{i2x2y,lp,mmd,pbrt}.txt \&\&

cmp resultados/i2x2y.txt resultados/pbrt.txt \&\&

echo "i2x2y e pbrt sao identicos."

```



\### Resultado



Cada método gerou 64 pontos distintos, com as duas coordenadas entre 0

inclusive e 1 exclusivo. As saídas de `i2x2y` e `pbrt` foram idênticas,

inclusive na ordem dos pontos. Os quatro arquivos estão na pasta `resultados`.



Na execução no Windows, o programa exibiu uma taxa de geração negativa.

Esse valor foi descartado como medida de desempenho.

