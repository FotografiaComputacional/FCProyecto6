%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Id completando el script y las funciones 
% que se piden pedidas en los diferentes apartados
clear; clc;
load puntos_distintivos

% Apartado 0) probar get_proy y error_T
disp("--------Apartado 0) probar get_proy y error_T--------")
xy1 = [ 200 300 400 500 600; 375 198 125 375 200];
xy2=[ 87 200 319 500 605; 445 186 69 375 100];

fprintf('Probando get_proy con puntos de prueba...\n');
T = get_proy(xy1, xy2);
fprintf('Matriz T obtenida:\n');
show_matriz(T);

% disp(aplica_T(xy1, T));

% Errores
errores = error_T(xy1, T, xy2);
fprintf('\nErrores para cada punto:\n');
for i = 1:length(errores)
    fprintf('Punto %d: %.4f\n', i, errores(i));
end

fprintf('\nError m�ximo: %.4f\n', max(errores));

% Apartado 1) Emparejamiento de puntos
%             Filtrado de parejas usando RANSAC
%             Estimaci�n de las transformaciones T

%puntos = S(5).xy;
%num_puntos = size(puntos, 2);
%img5 = imread('pano_05.jpg')
%figure;
%imshow(img5);
%hold on;
%plot(puntos(1,:), puntos(2,:), 'r.', 'MarkerSize', 10);
%title(sprintf('Imagen %d con %d puntos distintivos', 5, num_puntos));

fprintf('Emparejamientos con la imagen n� 5:\n');
for i = 1:length(S)
    if i ~= 5
        [xy1, xy2] = emparejar(S(i), S(5));
        
        num_emparejamientos = size(xy1, 2);
        fprintf('(%d,5) -> %d posibles parejas\n', i, num_emparejamientos);
    end
end


T = zeros(3, 3, 12);
T(:, :, 5) = eye(3);

fprintf('Emparejamientos con la imagen n� 5:\n');
for i = 1:length(S)
    if i ~= 5
        [xy1, xy2] = emparejar(S(i), S(5));
        
        num_emparejamientos = size(xy1, 2);
        fprintf('(%d, 5) -> %d posibles parejas', i, num_emparejamientos);
        
        [xy1_filtered, xy2_filtered] = ransac(xy1, xy2);
        
        num_aceptados = size(xy1_filtered, 2);
        fprintf(' -> %d aceptadas\n', num_aceptados);
        
        if num_aceptados >= 4
            T(:, :, i) = get_proy(xy1_filtered, xy2_filtered);
        end
    end
end

fprintf('\nMatriz de transformaci�n para la imagen 1:\n');
show_matriz(T(:, :, 1));

fprintf('\nMatriz de transformaci�n para la imagen 12:\n');
show_matriz(T(:, :, 12));

% Guardar el array T
save('T_data.mat', 'T');
fprintf('\nMatrices de transformaci�n guardadas en T_data.mat\n');



% Apartado 2) Montaje del mosaico 
%             Eliminaci�n de franjas negras

load('T_data.mat');

img_ref = im2double(imread('pano_01.jpg'));
[N, M, ~] = size(img_ref);

Umin = 1;
Umax = M;
Vmin = 1;
Vmax = N;

[~, ~, NF] = size(T);

for k = 1:NF
    fich = sprintf('pano_%02d.jpg', k);
    im = im2double(imread(fich));

    [~, RU, RV] = warp_img(im, T(:,:,k));
    
    Umin = min(Umin, RU(1));
    Umax = max(Umax, RU(2));
    Vmin = min(Vmin, RV(1));
    Vmax = max(Vmax, RV(2));

end

% Mostrar el rango obtenido para las coordenadas u y v
fprintf('Rango de coordenadas del mosaico: U=[%.1f, %.1f], V=[%.1f, %.1f]\n', Umin, Umax, Vmin, Vmax);

dU = (1 - Umin);
dV = (1 - Vmin);

anchoMosaico = Umax + dU;
altoMosaico = Vmax + dV;

fprintf('Dimensiones del mosaico final: %.1f x %.1f\n', anchoMosaico, altoMosaico);
msaico = zeros(ceil(altoMosaico), ceil(anchoMosaico), 3);

for k = 1:NF
    fich = sprintf('pano_%02d.jpg', k);
    im = im2double(imread(fich));
    
    [im_warped, RU, RV] = warp_img(im, T(:, :, k));
    
    col_ini = ceil(RU(1) + dU);
    row_ini = ceil(RV(1) + dV);
    col_fin = col_ini + size(im_warped, 2) - 1;
    row_fin = row_ini + size(im_warped, 1) - 1;

    [alto_m, ancho_m, ~] = size(msaico);

    % recortamos si se sale porque salia error 
    if row_fin > alto_m
        im_warped = im_warped(1:(alto_m - row_ini + 1), :, :);
        row_fin = alto_m;
    end
    if col_fin > ancho_m
        im_warped = im_warped(:, 1:(ancho_m - col_ini + 1), :);
        col_fin = ancho_m;
    end

    Q = msaico(row_ini:row_fin, col_ini:col_fin, :);
    for c = 1:3
        canal = im_warped(:, :, c);
        Qc = Q(:, :, c);
        Qc(~isnan(canal)) = canal(~isnan(canal));
        Q(:, :, c) = Qc;
    end
    msaico(row_ini:row_fin, col_ini:col_fin, :) = Q;
end


figure; 
imshow(msaico);
title('Mosaico apartado 2');


% Apartado 3) Ajuste GLOBAL del mosaico
%             Gr�ficas del error global antes/despues optimizacion

%P = find_solapamientos(S);
%save('P_data.mat', 'P');

load('P_data.mat');
N = sum([P.np]);
fprintf("Número total de puntos emparejados: %d\n", N);
disp(P(3))

ancla = 5;
param0 = T2param(T, ancla);
err = error_global(param0, P, ancla);
figure;
plot(err);
title('Error global antes de la optimización');
xlabel('Índice');
ylabel('Error (coordenadas u y v)');

s = std(err);
fprintf('Desviación estándar del error global (antes de optimización): σ = %.4f\n', s);

f_min=@(param)error_global(param,P,ancla);
param_def=lsqnonlin(f_min,param0);
err_opt = error_global(param_def, P, ancla);
figure;
plot(err_opt);
title('Error global después de la optimización');
xlabel('Índice');
ylabel('Error (coordenadas u y v)');
s = std(err_opt);
fprintf('Desviación estándar del error global (después de optimización): σ = %.4f\n', s);

T_opt = param2T(param_def, ancla);
fprintf('\nMatriz de transformación optimizada T(:,:,1):\n');
show_matriz(T_opt(:,:,1));

fprintf('\nMatriz de transformación optimizada T(:,:,12):\n');
show_matriz(T_opt(:,:,12));

%Cambio de imagen de referencia
NF = 12;         % Número de imágenes
ancla = 2   ;       % Imagen de referencia

T_malas = repmat(eye(3), [1, 1, NF]);
param_malo = T2param(T_malas, ancla);
err_ini = error_global(param_malo, P, ancla);

fprintf('---- Cambio de imagen de referencia -----:\n');

figure;
plot(err_ini);
title('Error global con hipótesis inicial incorrecta');
xlabel('Índice');
ylabel('Error');
fprintf('Desviación estándar del error inicial: σ = %.4f\n', std(err_ini));

f_malo = @(param) error_global(param, P, ancla);
param_corr = lsqnonlin(f_malo, param_malo);
err_corr = error_global(param_corr, P, ancla);

figure;
plot(err_corr);
title('Error global después de optimización (partiendo de identidad)');
xlabel('Índice');
ylabel('Error');

fprintf('Desviación estándar del error después de optimización: σ = %.4f\n', std(err_corr));

T_final = param2T(param_corr, ancla);

fprintf('\nMatriz de transformación T(:,:,1) con ancla en imagen 2:\n');
show_matriz(T_final(:, :, 1));

% Mosaico con imagen de referencia 2 - Igual que antes pero con T_final

%%%%%%%%  FIN del SCRIPT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_matriz(A)
 fprintf('%7.4f %7.4f %9.2f\n',A'); 
 fprintf('--------\n');
end

%%%%%%%%%%%% Funciones auxiliares a completar %%%%%%%%%%%%%%
%%%%%%% en los distintos apartados de la practica  %%%%%%%%%

% FUNCIONES Apartado 0)  Funciones get_proy(), aplica_T y error_T 

function T = get_proy(xy, uv)
    Np = size(xy, 2);
    
    x = xy(1, :)';
    y = xy(2, :)';
    u = uv(1, :)';
    v = uv(2, :)';

    A = zeros(2*Np, 8);
    b = zeros(2*Np, 1);
    
    for i = 1:Np
        % Ecuaciones para u
        A(2*i-1, 1) = x(i);
        A(2*i-1, 2) = y(i);
        A(2*i-1, 3) = 1;
        A(2*i-1, 7) = -u(i)*x(i);
        A(2*i-1, 8) = -u(i)*y(i);
        b(2*i-1) = u(i);
        
        % Ecuaciones para v
        A(2*i, 4) = x(i);
        A(2*i, 5) = y(i);
        A(2*i, 6) = 1;
        A(2*i, 7) = -v(i)*x(i);
        A(2*i, 8) = -v(i)*y(i);
        b(2*i) = v(i);
    end
    
    params = A \ b;
    
    % Parametros
    a = params(1);
    b = params(2);
    c = params(3);
    d = params(4);
    e = params(5);
    f = params(6);
    g = params(7);
    h = params(8);
    
    T = [a b c; d e f; g h 1];
end

function [uv] = aplica_T(xy, T)
    % Paso 1
    xy_homogeneas = [xy; ones(1, size(xy, 2))];
    % Paso 2
    resultado = T * xy_homogeneas;
    % Paso 3
    uv = resultado(1:2, :) ./ repmat(resultado(3, :), 2, 1);
end

function err = error_T(xy1, T, xy2)
    xy2_pred = aplica_T(xy1, T);
    
    diff = xy2_pred - xy2;
    err = sqrt(sum(diff.^2, 1));
end

% Apartado 1)  

function [k,d]=nearest(id,ID)
    DIF = ID - id;
    distancias = sqrt(sum(DIF.^2, 1));
    [d, k] = min(distancias);
end

function [xy1,xy2]=emparejar(s1,s2)
    ID1 = s1.id;
    ID2 = s2.id;
    
    N1 = size(ID1, 2);
    
    pareja_de = zeros(1, N1);
    
    % Bucle en k barriendo los N1 puntos
    for k = 1:N1
        [j, d] = nearest(ID1(:,k), ID2);
        [k_nearest, d_nearest] = nearest(ID2(:,j), ID1);
        
        if k_nearest == k && d < 0.5
            pareja_de(k) = j;
        end
    end
    
    idx1 = find(pareja_de > 0);
    idx2 = pareja_de(idx1);
    xy1 = s1.xy(:, idx1);
    xy2 = s2.xy(:, idx2);
end


function [XY1, XY2] = ransac(XY1, XY2)
    % 1) Determinar N
    N = size(XY1, 2);
    
    % 2) Inicializar un vector ok
    ok = [];
    
    % 3) Hacer un bucle de 2000 iteraciones
    for iter = 1:2000
        % Elegir 4 de los N puntos al azar
        idx = randperm(N, 4);
        
        T = get_proy(XY1(:, idx), XY2(:, idx));
        errores = error_T(XY1, T, XY2);
        cercanos = find(errores < 1);
            
        if length(cercanos) > length(ok)
            ok = cercanos;
        end
    end
    
    % 4) �ndices de los puntos
    XY1 = XY1(:, ok);
    XY2 = XY2(:, ok);
end


% FUNCIONES del Apartado 2)

function [im2, RU, RV] = warp_img(im, T)
 % 1) Convertir la imagen a double (valores entre 0 y 1) y hallar su alto N y ancho M
 im = im2double(im);
 [N, M, ~] = size(im);
 
 % 2) Construir los vectores RU y RV con el rango de coordenadas de la imagen en el espacio destino
 x = [1, M, M, 1];
 y = [1, 1, N, N];
 xy = [x; y];
 
 uv = aplica_T(xy, T);
 
 u_min = floor(min(uv(1,:)));
 u_max = ceil(max(uv(1,:)));
 v_min = floor(min(uv(2,:)));
 v_max = ceil(max(uv(2,:)));
 
 RU = u_min:u_max;
 RV = v_min:v_max;
 
 % 3) Reservar memoria para la imagen de salida im2
 N_new = length(RV);
 M_new = length(RU);
 im2 = zeros(N_new, M_new, 3);

 X = zeros(N_new, M_new);
 Y = zeros(N_new, M_new);
 
 % 4) Calcular T^-1 la matriz inversa de T (para el warping "inverso")
 T_inv = inv(T);
 
 % 5) Hacer un bucle para barrer todas las filas de la imagen destino
 for k = 1:N_new
     % a) Crear matriz uv con las coordenadas u y v de esa fila
     u = RU;
     v = RV(k) * ones(size(RU));
     uv = [u; v];
     % b) Obtener las correspondientes coordenadas xy en la imagen de partida
     xy = aplica_T(uv, T_inv);
     % c) Guardar las coordenadas
     X(k,:) = xy(1,:);
     Y(k,:) = xy(2,:);
 end
 
 % 6) Usar interp2 para interpolar los valores de la imagen en las coordenadas (X,Y)
 for c = 1:3
     im2(:,:,c) = interp2(1:M, 1:N, im(:,:,c), X, Y, 'linear');
 end
end


% FUNCIONES de Apartado 3)

function param = T2param(T, ancla)
    NF = size(T, 3);         % Número de ima
    param = [];              

    for k = 1:NF
        if k ~= ancla
            M = T(:,:,k);
            p = reshape(M', 1, 9);  %' para ir filaa fila
            param = [param, p(1:8)];
        end
    end
end

function T = param2T(param, ancla)
    NF = length(param)/8 + 1;  % el número de im
    T = zeros(3, 3, NF);       

    idx = 1;
    for k = 1:NF
        if k == ancla
            T(:,:,k) = eye(3); 
        else
            p = param(idx:idx+7);
            M = [p(1) p(2) p(3);
                 p(4) p(5) p(6);
                 p(7) p(8) 1];
            T(:,:,k) = M;
            idx = idx + 8;
        end
    end
end


function P = find_solapamientos(S)
    NF = length(S);
    k = 1;                  

    for i = 1:NF
        for j = i+1:NF
            [xy1, xy2] = emparejar(S(i), S(j));

            if size(xy1, 2) >= 10
                % Filtrado con RANSAC
                [xy1_f, xy2_f] = ransac(xy1, xy2);
                np = size(xy1_f, 2);

                if np >= 10
                    P(k) = struct( ...
                        'i', i, ...
                        'j', j, ...
                        'np', np, ...
                        'xy1', xy1_f, ...
                        'xy2', xy2_f);
                    k = k + 1;
                end
            end
        end
    end
end

function err = error_global(param, P, ancla)
    T = param2T(param, ancla);
    N = sum([P.np]);
    dif = zeros(2, N);

    rg = 0;             

    for k = 1:length(P)
        np = P(k).np;                    
        
        rg = rg(end) + (1:np); 

        i = P(k).i;
        j = P(k).j;

        uv1 = aplica_T(P(k).xy1, T(:, :, i));
        uv2 = aplica_T(P(k).xy2, T(:, :, j));

        dif(:, rg) = uv2 - uv1;
    end

    err = reshape(dif, 1, []); 
end










