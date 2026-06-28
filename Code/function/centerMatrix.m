function B = centerMatrix(A, M, N)
    % 将矩阵A居中放置在MxN的零矩阵中
    B = zeros(M, N);
    [a_rows, a_cols] = size(A);

    % 确保目标矩阵足够大
    if a_rows > M || a_cols > N
        error('目标矩阵太小，无法容纳原矩阵');
    end

    % 计算起始位置
    r_start = floor((M - a_rows)/2) + 1;
    c_start = floor((N - a_cols)/2) + 1;

    % 放置矩阵
    B(r_start:r_start+a_rows-1, c_start:c_start+a_cols-1) = A;
end