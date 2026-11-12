# Análise de amostragem Sobol 2D

Projeto de pesquisa de Computação Gráfica — IFNMG, 2026/2.

Investigamos se LP e MMD reduzem o erro da renderização em relação à Sobol 2D convencional, mantendo as condições do experimento iguais. Ainda não há resultado quantitativo de qualidade.

## Origem dos códigos

- `benchmarking.cpp` e `sobolmatrices.h`: suplemento de Abdalla G. M. Ahmed, *An Implementation Algorithm of 2D Sobol Sequence: Fast, Elegant, and Compact*, EGSR 2024. DOI: https://doi.org/10.2312/sr.20241147
- Suplemento oficial, arquivo `paper1031_supp.zip`: https://diglib.eg.org/items/57f2cdeb-69d9-434e-8cf8-37b63e7e69d9
- `renderizador/smallpt.cpp`: smallpt original de Kevin Beason, https://www.kevinbeason.com/smallpt/ . Aviso MIT em `renderizador/LICENCA_smallpt.txt`.
- `renderizador/amostragem.h`: funções J e P extraídas do suplemento de Ahmed, mais as funções de sementes e deslocamento digital usadas na integração.
- `renderizador/smallpt_pesquisa.cpp`: adaptação do smallpt para consumir os pontos Sobol nas duas coordenadas da amostragem da câmera.

O OpenQMC (https://github.com/AcademySoftwareFoundation/openqmc) foi consultado como referência técnica. Não é dependência desta implementação.

## Ambiente de execução

Os comandos abaixo devem ser executados no terminal **MSYS2 MSYS**, na pasta do repositório. O caminho da pasta varia conforme o computador. Caso o compilador ainda não esteja instalado nesse ambiente:

```bash
pacman -S --needed gcc
```

Os programas usam funções POSIX, como `drand48` e `erand48`. A execução no MSYS2 MSYS foi escolhida após uma dificuldade com `drand48` no UCRT64. Os executáveis não são versionados; devem ser compilados localmente.

## Teste de 64 pontos

Compile o programa do suplemento e gere os arquivos:

```bash
g++ -std=gnu++17 -O2 -include ctime benchmarking.cpp -o benchmarking.exe
mkdir -p resultados
for metodo in i2x2y lp mmd pbrt; do
  ./benchmarking.exe -a "$metodo" -m 6 "resultados/$metodo.txt" || break
done
```

`-m 6` pede 2 elevado a 6, ou 64 pontos. A primeira linha informa a quantidade e as demais contêm pares de coordenadas.

Para validar a quantidade, a unicidade e o intervalo das coordenadas:

```bash
awk '
  FNR == 1 { if ($1 != 64) erro = 1; next }
  { total[FILENAME]++ }
  NF != 2 || $1 < 0 || $1 >= 1 || $2 < 0 || $2 >= 1 { erro = 1 }
  !vistos[FILENAME, $1 + 0, $2 + 0]++ { distintos[FILENAME]++ }
  END {
    for (i = 1; i < ARGC; i++) {
      arquivo = ARGV[i]
      printf "%s: %d pontos, %d distintos\n", arquivo, total[arquivo], distintos[arquivo]
      if (total[arquivo] != 64 || distintos[arquivo] != 64) erro = 1
    }
    print erro ? "Falha na validacao." : "Todos os pontos estao no intervalo esperado."
    exit erro
  }
' resultados/{i2x2y,lp,mmd,pbrt}.txt &&
cmp resultados/i2x2y.txt resultados/pbrt.txt &&
echo "i2x2y e pbrt sao identicos."
```

Cada método gerou 64 pontos distintos em `[0,1)²`. As saídas de `i2x2y` e `pbrt` foram idênticas, inclusive na ordem. Os arquivos estão em `resultados/`. A taxa de geração negativa exibida nessa execução foi descartada como medida de desempenho.

## Renderização original

```bash
cd renderizador
g++ -std=gnu++98 -O2 smallpt.cpp -o smallpt.exe
./smallpt.exe 4
```

O original gera `image.ppm` em 1024 × 768 com quatro amostras por pixel. Essa execução serviu para verificar o funcionamento do renderizador antes da adaptação. Sua imagem não é uma comparação controlada com o teste abaixo, pois resolução, orçamento de amostras e fluxo de sementes são diferentes.

## Primeira renderização com Sobol

Execute dentro de `renderizador`:

```bash
g++ -std=gnu++17 -O2 smallpt_pesquisa.cpp -o smallpt_pesquisa.exe
./smallpt_pesquisa.exe sobol 16 sobol_16.ppm 320 240 1
```

A ordem dos argumentos é método, amostras por pixel, arquivo de saída, largura, altura e semente. A quantidade de amostras por pixel deve ser uma potência de dois a partir de quatro.

O resultado executado no Windows é `renderizador/sobol_16.ppm`: 320 × 240, 16 amostras por pixel e semente 1. A imagem foi aberta e conferida. Esse é um teste de funcionamento da integração, sem conclusão sobre melhoria de qualidade.

### O que a integração modifica

O smallpt divide cada pixel em quatro subpixels. As 16 amostras por pixel correspondem a quatro pontos Sobol por subpixel. Suas duas coordenadas passam pela mesma transformação de filtro tent do original e escolhem a posição dos raios da câmera. Por causa desse filtro, o suporte da amostragem pode ultrapassar os limites de um subpixel.

Um deslocamento digital por XOR, fixo para cada subpixel e semente, evita repetir o mesmo padrão em toda a imagem. Esse deslocamento não é embaralhamento de Owen. Cena, câmera, materiais e função de transporte de luz foram mantidos. As direções de reflexão e as demais decisões do transporte continuam usando `erand48`.

A adaptação separa os estados pseudoaleatórios da câmera e do transporte e define sementes por amostra. Portanto, não reproduz bit a bit o fluxo de sementes por linha do original.

Há também um controle auxiliar aleatório no mesmo programa:

```bash
./smallpt_pesquisa.exe random 16 random_16.ppm 320 240 1
```

Esse comando está disponível, mas sua execução pela equipe ainda não foi registrada. A baseline da pesquisa é Sobol, e as variantes a integrar são LP e MMD.

## Checkpoint 2

A baseline foi reproduzida, o teste de 64 pontos foi registrado e a primeira
renderização com Sobol foi concluída. O [esqueleto do artigo](artigo/artigo.pdf)
foi compilado no Overleaf. O [arquivo LaTeX](artigo/artigo.tex) e as
[referências](artigo/referencias.bib) estão no repositório.

LP e MMD ainda não foram integrados ao renderizador.

## Após o checkpoint

Concluir LP e MMD, respeitando as condições e os tamanhos dos conjuntos do artigo. Comparar os três métodos com mesma cena, câmera, resolução, filtro, amostras por pixel e política de sementes, usando repetições. O escopo inicial altera somente o par de dimensões da câmera.

Antes dos experimentos quantitativos, preparar saída de radiância linear sem o recorte por subpixel do smallpt e sem correção gama. O PPM atual mantém o recorte e a transformação de exibição do original e serve para inspeção visual; não será usado diretamente para RMSE linear. Gerar e validar uma referência de alta amostragem e então medir erro e convergência. Não há evidência atual de superioridade de LP ou MMD.
