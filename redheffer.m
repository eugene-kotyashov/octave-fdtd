function res = redheffer(sa11,sa12,sa21,sa22,sb11,sb12,sb21,sb22)
	D = sa12*(eye(2)-sb11*sa22)^(-1);
	F = sb21*(eye(2)-sa22*sb11)^(-1);
	s11 = sa11 + D*sb11*sa21;
	s12 = D*sb12;
	s21 = F*sa21;
	s22 = sb22 + F*sa22*sb12;
	res = [s11;s12;s21;s22];
endfunction
