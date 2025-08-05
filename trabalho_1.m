%Código para o Trabalho 1 da disciplina SME 0205, Métodos do Cálculo Numérico I

% Item 0: Código para receber os arquivos manh.el e manh.xy e plotar o grafo de ruas da ilha de Manhattan, NY.
edgeFile = 'manh.el';
edges = dlmread(edgeFile);
coordFile = 'manh.xy';
coords = dlmread(coordFile);
numVertices = max(max(edges));
adjMatrix = zeros(numVertices+1, numVertices+1);
for i = 1:size(edges, 1)
 vertex1 = edges(i, 1) + 1;
 vertex2 = edges(i, 2) + 1;
 adjMatrix(vertex1, vertex2) = 1;
 adjMatrix(vertex2, vertex1) = 1;
end
figure;
gplot(adjMatrix, coords, '-o');
title('Grafo das Ruas');
xlabel('Longitude');
ylabel('Latitude');

% Item 1: Código para selecionar os vértices aleatoriamente e atribuir valores neles
graph = digraph(adjMatrix);
connectedComponents = conncomp(graph);
largestComponentIndex = mode(connectedComponents);
largestComponentIndices = find(connectedComponents == largestComponentIndex);
k = 800;
selectedVertices = randsample(largestComponentIndices, k);
selectedValues = randi([0, 800], k, 1);
figure;
gplot(adjMatrix, coords, '-o');
hold on;
plot(coords(selectedVertices, 1), coords(selectedVertices, 2), 'ro', 'MarkerSize', 10);
text(coords(selectedVertices, 1), coords(selectedVertices, 2), num2str(selectedValues), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off;
title('Grafo com Vértices Selecionados');
xlabel('Longitude');
ylabel('Latitude');
% Item 2: Código para construir a matriz Laplaciana L do grafo das ruas
degreeMatrix = diag(sum(adjMatrix, 2));
L = degreeMatrix - adjMatrix;
figure;
imagesc(L);
colorbar;
title('Matriz Laplaciana');
xlabel('Vértices');
ylabel('Vértices');
% Item 3: Código para construir a matriz de penalidades P
alpha = 1.0e7;
P = zeros(numVertices+1, numVertices+1);
P(selectedVertices, selectedVertices) = alpha * eye(k);
figure;
imagesc(P);
colormap('jet');
colorbar;
title('Matriz de Penalidades');
% Item 4: Código para construir o vetor b
b = zeros(numVertices+1, 1);
for i = 1:k
 vertexIndex = selectedVertices(i);
 b(vertexIndex) = selectedValues(i);
end
figure;
stem(b);
title('Vetor b');
% Item 5: Código para resolver o sistema: (matriz Laplaciana + matriz de penalidades) * x = matriz de penalidades * vetor b
% Método 1: Decomposição LU
tic;
A_lu = L + P;
x_lu = A_lu\b;
time_lu = toc;
% Método 2: Jacobi
tic;
D = diag(diag(L));
R = L - D;
M_jacobi = -D\R;
N_jacobi = D\b;
x_jacobi = zeros(numVertices+1, 1);
maxIterations = 1000;
tolerance = 1e-6;
for iter = 1:maxIterations
   x_jacobi_new = M_jacobi * x_jacobi + N_jacobi;
   if norm(x_jacobi_new - x_jacobi) < tolerance
       break;
   end
   x_jacobi = x_jacobi_new;
end
time_jacobi = toc;
% Método 3: Gauss-Seidel
tic;
D = diag(diag(L));
L_gs = tril(L);
U_gs = triu(L);
M_gs = -(D + L_gs)\U_gs;
N_gs = (D + L_gs)\b;
x_gs = zeros(numVertices+1, 1);
for iter = 1:maxIterations
   x_gs_new = M_gs * x_gs + N_gs;
   if norm(x_gs_new - x_gs) < tolerance
       break;
   end
   x_gs = x_gs_new;
end
time_gs = toc;
% Método 4: Gradientes Conjugados
tic;
A_gc = L + P;
x_gc = pcg(A_gc, b, tolerance, maxIterations);
time_gc = toc;
% Comparação dos tempos de solução
fprintf('Tempo de solução (LU): %.6f segundos\n', time_lu);
fprintf('Tempo de solução (Jacobi): %.6f segundos\n', time_jacobi);
fprintf('Tempo de solução (Gauss-Seidel): %.6f segundos\n', time_gs);
fprintf('Tempo de solução (Gradientes Conjugados): %.6f segundos\n', time_gc);
% Visualização dos resultados
figure;
subplot(2, 3, 1);
stem(x_lu);
title('Resultado (LU)');
subplot(2, 3, 2);
stem(x_jacobi);
title('Resultado (Jacobi)');
subplot(2, 3, 3);
stem(x_gs);
title('Resultado (Gauss-Seidel)');
subplot(2, 3, 4);
stem(x_gc);
title('Resultado (Gradientes Conjugados)');
% Imagens das soluções por cada método
% Método 1: LU
figure;
stem(x_lu);
title('Resultado (LU)');
xlabel('Vértices');
ylabel('Valor');
saveas(gcf, 'resultado_lu.png');
% Método 2: Jacobi
figure;
stem(x_jacobi);
title('Resultado (Jacobi)');
xlabel('Vértices');
ylabel('Valor');
saveas(gcf, 'resultado_jacobi.png');
% Método 3: Gauss-Seidel
figure;
stem(x_gs);
title('Resultado (Gauss-Seidel)');
xlabel('Vértices');
ylabel('Valor');
saveas(gcf, 'resultado_gauss_seidel.png');
% Método 4: Gradientes Conjugados
figure;
stem(x_gc);
title('Resultado (Gradientes Conjugados)');
xlabel('Vértices');
ylabel('Valor');
saveas(gcf, 'resultado_gradientes_conjugados.png');
% Tabela comparando os tempos de solução
methods = {'LU', 'Jacobi', 'Gauss-Seidel', 'Gradientes Conjugados'};
times = [time_lu, time_jacobi, time_gs, time_gc];
tableData = table(methods', times', 'VariableNames', {'Método', 'Tempo de Solução'});
disp(tableData);
