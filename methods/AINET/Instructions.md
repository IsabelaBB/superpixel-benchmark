Para que funcione conforme o driver NVIDIA instalado, execute:
1) pipenv install # instala o ambiente
2) nvidia-smi # copie a versão cuda encontrada
2) pipenv shell 
3) Execute o comando indicado para sua versão em https://pytorch.org/get-started/previous-versions. Exemplo para versão cuda 10.1: pip install torch==1.8.1+cu101 torchvision==0.9.1+cu101 torchaudio==0.8.1 -f https://download.pytorch.org/whl/cu101/torch_stable.html
4) Teste se a versão instalada reconhece o driver cuda: $ python3 -c "import torch; print(torch.cuda.is_available())"