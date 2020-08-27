# -*- coding: utf-8 -*-
"""
Created on Sat May  4 12:17:06 2019

@author: Marcelo Rocha
"""
import numpy as np
import math

#OBS: A e B sao elementos de W usados como referencia
def calc_c_s(A,B):
    a=b=0
    if A<0:  #ver o modulo
        a=A*(-1)
    if B<0:
        b=B*(-1)
    if a>b:
        if A !=0 :
            t=-B/A
            c=1/math.sqrt(1+t**2)
            s=c*t
        if A==0:
            c=s=0
    if a<=b:
        if  B !=0:
            t=-A/B
            s=1/math.sqrt(1+t**2)
            c=s*t
        if B==0:
            c=s=0
    return c,s
    
#OBS:MATRIZ W com valores FLUTUANTES
#OBS: i é a coluna e j é a linha do elemento a ser modificado !
def rotgivens(W,m,i,j,c,s):   #vai ser chamada na main para alterar elemento por elemento atraves de um for que varre as colunas
    R=[[]]        #Cria a matriz que vira a ser a triangular superior
    R=np.copy(W)  #Copia a matriz W em R,inicialmente 
   # i=i-1  #ajuste ao indice
   # j=j-1  #ajuste ao indic

    for p in range(m):  #varre as colunas de W
        R[i][p]=c*W[i][p]-s*W[j][p]
        R[j][p]=s*W[i][p]+c*W[j][p]
    return R

def transfb(R,W,b): #função que altera o vetor b para btransf,recebe como parametros a matriz triangular R,a matriz W e o vetor b a ser modificado
    Q=[[]]          #matriz ortogonal de transformação
    x=[]           #solução de Rx=btransf
    R_inversa=np.linalg.pinv(R)
    Q=np.dot(W,R_inversa)
    Q_inversa=np.linalg.pinv(Q)
    btransf=np.dot(Q_inversa,b)
    x=np.dot(R_inversa,btransf)
    return x
                            
                           #Varios sistemas simultaneos
def sistsimultaneos(A,W,n,m,p): #Atraves dos parametros m,n,p de W e A,transforma a matriz A e W pela rotação de givens e,alem  disso ,obtem H,solução da minimização do erro
    H=[[]]    #h é a matriz que resolve o sistema A=W*H,logo,H=W_inversa*A. Posteriormente sera armazenada a solução 
    W_inversa=np.linalg.pinv(W)
    #print('n:',n,'m:',m,'p:',p)
    H=np.dot(W_inversa,A)
    W_transf=W  # tanto faz copiar ou pegar a mesma
    A_transf=A
    somatoria=0 #valor que armazena a soma do produto de WH
  

    
    for k in range(p):   #varre coluna a coluna da matriz. Laço que transforma A e W em matrizes triangulares superiores
        for j in range(n-1,k,-1):  #varre a coluna,elemento a elemento até chegar na diagonal principal (j=k)
            i=j-1
            if W_transf[j][k] !=0 and A_transf[j][k] !=0:
                c,s=calc_c_s(W_transf[i][k],W_transf[j][k])     
                W_transf=rotgivens(W_transf,p,i,j,c,s)  #aplica rotação de givens à W
                #cos,sen=calc_c_s(A_transf[i][k],A_transf[j][k])   #####
                A_transf=rotgivens(A_transf,m,i,j,c,s) #aplica rotação de givens à A
    p=p-1 # ajuste pra p poder ser ajustado a indice,parece se ter uma comparação entre indices e diminuir a confusão     
    n=n-1 # ajuste de n pra indice para comparar valores entre indices em python
    for k in range(p,-1,-1): # vai de p a 0,define os elementos de H
        for j in range(m):
            for x in range(k,k+p): #laço que realiza a somatoria
                for y in range(k+1,p+1):
                    if k+p-1<=n and k+1<=p and y==x+1:
                       somatoria+=W[x][y]*H[y][j]
                       
            if W[k][k] !=0:           
                H[k][j]=(A[k][j]-somatoria)/W[k][k] #obs: pra cada elemento H[k][j],é realizado o laço da somatoria
            somatoria=0
    
    return H


def fatnnegativa(A,n,m,p):  #função que fatora por matrizes nao negativas
    W=np.floor(10*np.random.random((n,p)))  #inicializa W randomicamente com valores entre 0 a 10
    A_copia=np.copy(A)
    s=sj=0
   # A_transposta=A.T
    W_transposta=[[]]
    
    for i in range(100): #estabilização da norma do erro,realização de itmax=100 iterações
        #print('fatnegat:',i)
        
        
        for j in range(p):  #Normalização de W: varre as colunas de W pra normaliza-la
            for lin in range(n): #determinação de sj
                s+=(W[lin][j])**2
                if lin==n-1:   
                    sj=math.sqrt(s)  ########
            for lin in range(n): #normalização de cada coluna
                if sj !=0:
                    W[lin][j]=W[lin][j]/sj
            s=sj=0
            
        H=sistsimultaneos(A_copia,W,n,m,p)
        for col in range(m): # laço que redfine H,começa varrendo as colunas de H
            for lin in range(p):
                if H[lin][col]<0:
                   H[lin][col]=0
    
        A_transpcopi=np.copy(A)
        A_transposta=A_transpcopi.T ###
        H_transposta=H.T
        W_transposta=sistsimultaneos(A_transposta,H_transposta,m,n,p)
        W=W_transposta.T
        for col in range(p):  #redefine W
            for lin in range(n):
                if W[lin][col]<0:
                    W[lin][col]=0
    return W,H

def matrizA(nome_arquivo):  
    arquivo=open(nome_arquivo,'r')
    linha=arquivo.readlines()
    arquivo.close()
    matriz=[[]]
    linhamat=[]
    e=''
    for i in linha:
        #print(i)
        for j in i:
            #print(j)
    
            if j !=' ' and j !='\n':
                e+=j
            if j ==' ' or j=='\n' and e !='':
              linhamat.append(float(e)/255)

              e=''
            if j == '\n':
                matriz.append(linhamat)
                linhamat=[]
    
      #  linhamat.append(e)
      #  matriz.append(linhamat)
     
    del matriz[0]    # primeiro indice é vazio

    return matriz

def main():
    nome_imagens=['train_dig0.txt','train_dig1.txt','train_dig2.txt','train_dig3.txt','train_dig4.txt','train_dig5.txt','train_dig6.txt','train_dig7.txt','train_dig8.txt','train_dig9.txt'] # Cada imagem corresponde a um A(um digito original) 
    matrizestreinadas=[]  # Armazena as matrizes W, w0,w1...w9
    p=10  # decompor em p imagens.Não é dado nem obtido, é definido

    for i in range(len(nome_imagens)):  # para cada digito.Laço que realiza o treinamento 
        print('treinamento da imagem:',i)
        m=n=0
        A=matrizA(nome_imagens[i])
        for linha in A:
            n+=1
            m=0
            for coluna in linha:
                m+=1                              
        W,H=fatnnegativa(A,n,m,p)
        matrizestreinadas.append(W)

    
    
    
    ########################### Laço que classifica o indice de acertos

    A=matrizA('test_images.txt')
    for linha in A:
        n_teste=0
        for coluna in linha:
            n_teste+=1  #numero de colunas de A
    
    digitos=[]
    for j in range(n_teste):       # para cada coluna j de A
        vetordigitos=[]
        vetorerros=[]
        print('Analise da coluna ',j,' da matriz teste')
        for n_w in range(len(matrizestreinadas)):  #varredura nos W com objetivo final achar o com menor erro correspondente à imagem J de A
            H=sistsimultaneos(A,matrizestreinadas[n_w],784,n_teste,p)
            C=np.subtract(A,np.dot(matrizestreinadas[n_w],H)) #matriz erro
            somatoria=erroj=0
            for i in range(784):
                somatoria+=C[i][j]**2
                if i==783:
                    print('n_w:',n_w)
                    erroj=math.sqrt(somatoria)
                    if n_w==0:
                       # print('Entrou')
                        vetorerros.append(float(erroj))
                        vetordigitos.append(float(n_w))
                    elif n_w>=1:
                        #print(('n_w:',n_w))
                        #print('comprimento da lista de treinadas:',len(matrizestreinadas))
                        #print('Vetor erros:',vetorerros)
                        if erroj<vetorerros[n_w-1]:
                            vetorerros.append(float(erroj))
                            vetordigitos.append(float(n_w))
                        else:
                            vetorerros.append(float(vetorerros[n_w-1]))
                            vetordigitos.append(float(n_w))
        digitos.append(vetordigitos[len(vetordigitos)-1]) #lista com os digitos.O elemento j dessa lista corresponde ao digito da coluna j de A
        
    #####depois de obtido os digitos,serao realizados os testes de assertividade
    arquivo=open('test_index.txt')
    l=arquivo.readlines()
    arquivo.close()
    e=''
    valoresreais=[]  #digitos reais das colunas de A
    for lh in l:
        for elemento in lh:
            if elemento !='\n':
                e+=elemento
            if elemento=='\n':
                valoresreais.append(float(e))
                e=''
    contador=[0,0,0,0,0,0,0,0,0,0]
    total=[0,0,0,0,0,0,0,0,0,0]  #numero de vezes que cada elemento apareceu,armazenado em total[elemento]
    for indice in range(10):   # Laço que armazena na lista contador o numero de vezes que o elemento indice se repete e amarzena na posição contador[indice]
        for alcance in range(len(digitos)):
            if digitos[alcance]==valoresreais[alcance]==indice:
                contador[indice]+=1
    v=-1
    for valor in valoresreais:
        v+=1
        for variavel in range(10):
            if valor==variavel:
                total[variavel]+=1
            
    totalacertos=0
    for i in contador:
        totalacertos+=i
    quantidadetotal=len(valoresreais)
    print('Porcentagem de acertos:',totalacertos/quantidadetotal)
    
    for i in range(10):
        if total[i] !=0:
            t=(contador[i]/total[i])*100
        print('Quantidade de acertos do digito ',i,':',contador[i])
        print('Porcentagem de acertos do digito',i,':',t)
    
        

if __name__ == "__main__": main()

            
    

                
        

       