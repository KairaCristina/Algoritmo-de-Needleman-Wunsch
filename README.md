# Algoritmo-de-Needleman-Wunsch
Algoritmo para alinhamento de sequências de nucleotídeos ou proteínas

```
def needleman_wunsch(seq1, seq2, match=1, mismatch=0, gap=0):
    n = len(seq1)
    m = len(seq2)
```


# Inicialização da matriz de pontuação

```
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

```

# Inicialização da matriz de direção

```
    direction_matrix = [[0] * (m + 1) for _ in range(n + 1)]

```

# Inicialização da primeira coluna e linha

```
    for i in range(n + 1):
        score_matrix[i][0] = gap * i
        direction_matrix[i][0] = 'UP'
    for j in range(m + 1):
        score_matrix[0][j] = gap * j
        direction_matrix[0][j] = 'LEFT'

```

# Preenchimento das matrizes
    
```
for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch

            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap

            max_score = max(diagonal_score, up_score, left_score)

            score_matrix[i][j] = max_score

            if max_score == diagonal_score:
                direction_matrix[i][j] = 'DIAG'
            elif max_score == up_score:
                direction_matrix[i][j] = 'UP'
            else:
                direction_matrix[i][j] = 'LEFT'

```
# Traçar o caminho de volta

```
    align1 = ''
    align2 = ''
    i = n
    j = m
    while i > 0 or j > 0:
        if direction_matrix[i][j] == 'DIAG':
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif direction_matrix[i][j] == 'UP':
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, score_matrix[n][m]

```

# Substitua a 'seq1' e 'seq2' pelas sequencias que deseja alinhar:

```
seq1 = "ATCGTAC"
seq2 = "ATGTTAT"
alignment1, alignment2, score = needleman_wunsch(seq1, seq2)
print("Alinhamento 1:", alignment1)
print("Alinhamento 2:", alignment2)
print("Pontuação do alinhamento:", score)

```
