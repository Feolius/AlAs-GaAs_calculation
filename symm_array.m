%Вспомогательная функция для симметризации решения. Создает вторую
%половину массива как "отраженную" первую
function arr_out = symm_array(arr)
global N
arr_out = zeros(1, N);
if rem(N, 2) == 1
    arr_out(1:(fix(N/2) + 1)) = arr(1:(fix(N/2) + 1));
    for g = (fix(N/2) + 2):N
        arr_out(g) = arr(N + 1 - g);        
    end
else
    arr_out(1:N/2) = arr(1:N/2);
    for g = (N/2 + 1):N
        arr_out(g) = arr(N + 1 - g);        
    end
end
