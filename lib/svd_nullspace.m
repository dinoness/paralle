function N = svd_nullspace(A, tol)
%SVD_NULLSPACE 用SVD求零空间，避免不同MATLAB版本下 null(A) 阈值差异较大
% 输入：
%   A   : m x n 矩阵
%   tol : 判零阈值，缺省时自动按矩阵大小估计
% 输出：
%   N   : n x k 零空间基，列向量构成零空间基

    if nargin < 2 || isempty(tol)
        tol = max(size(A)) * eps(norm(A));
    end

    [~, S, V] = svd(A);
    s = diag(S);

    if isempty(s)
        N = eye(size(A, 2));
        return;
    end

    r = sum(s > tol);
    N = V(:, r+1:end);
end
