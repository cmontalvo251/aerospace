function H = Cholesky(S)

[r,c] = size(S);

H = 0.*S;

H(1,1) = sqrt(S(1,1));

if r > 1
  H(2,1) = S(1,2)/H(1,1);
  H(2,2) = sqrt(S(2,2)-H(2,1)*conj(H(2,1)));
end
if r > 2
  H(3,1) = S(3,1)/H(1,1);
  H(3,2) = 0;
  if abs(H(2,2)) > 1e-8
    H(3,2) = (S(3,2)-H(3,1)*conj(H(2,1)))/H(2,2);
  end
  H(3,3) = sqrt(S(3,3)-H(3,1)*conj(H(3,1))-H(3,2)*conj(H(3,2)));
end
