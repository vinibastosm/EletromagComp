from PIL import Image
import os
import re

# Configurações
diretorio_imagens = "imagens_heatmap"  # Altere para o diretório que contém as imagens
nome_arquivo_gif = "erro.gif"
duracao = 100  # Duração de cada frame em milissegundos
loop = 0       # Loop infinito (0), ou número de vezes para repetir o GIF

# Função para ordenar por números dentro dos nomes dos arquivos
def ordenar_numerico(nome):
    return [int(texto) if texto.isdigit() else texto for texto in re.split(r'(\d+)', nome)]

# Coleta e ordena todas as imagens na pasta
imagens = [Image.open(os.path.join(diretorio_imagens, arquivo))
           for arquivo in sorted(os.listdir(diretorio_imagens), key=ordenar_numerico) if arquivo.endswith(".png")]

# Cria e salva o GIF
if imagens:
    imagens[0].save(
        nome_arquivo_gif,
        save_all=True,
        append_images=imagens[1:],  # Exclui a primeira imagem da lista adicional
        duration=duracao,
        loop=loop
    )
    print(f"GIF criado com sucesso e salvo como {nome_arquivo_gif}")
else:
    print("Nenhuma imagem encontrada no diretório.")
