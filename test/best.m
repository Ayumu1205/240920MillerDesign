% 目的関数の定義
function fun = objective(x)
  % a,b,thetaの設定  fminconm関数では、x(1), x(2), x(3)が変数として扱われるため、a,b,thetaをx(1), x(2), x(3)に代入
  a = x(1);
  b = x(2);
  theta = x(3);

  % 焦点座標を求める関数
  [f1x, f1y, f2x, f2y, f] = solveFocus(a, b);

  % M3座標を求める関数
  [M3x, M3y, M3] = solveM3(f1x, f2x, theta);

  % Wの座標 楕円と直線Lの交点を求める関数
  [crossW1, crossW2] = solveM1Intersections(a, b, f1x, f1y, theta);

  % 2つの交点のうち、M3に近い方を選択し、M1とする関数
  [W1] = judgePointFunction(crossW1, crossW2, M3);

  % W1とM3の距離
  W1M3 = sqrt((W1(1) - M3(1))^2 + (W1(2) - M3(2))^2);

  % 目的関数の定義  W1M3の最小化を試みる
  fun = W1M3;
    
end

% 制約関数の定義 (x > 0 の制約)
function [c, ceq] = constraint(x)   
	a = x(1);
	b = x(2);
	theta = x(3);

	% 焦点座標を求める関数
	[f1x, f1y, f2x, f2y, f] = solveFocus(a, b);

	% M3座標を求める関数
	[M3x, M3y, M3] = solveM3(f1x, f2x, theta);

	% Wの座標 楕円と直線Lの交点を求める関数
	[crossW1, crossW2] = solveM1Intersections(a, b, f1x, f1y, theta);

	% 2つの交点のうち、M3に近い方を選択し、M1とする関数
	[W1] = judgePointFunction(crossW1, crossW2, M3); 

	% Hの座標 楕円と直線Nの交点を求める関数
	[crossH1, crossH2] = solveH1Intersections(a, b, M3x, M3y, theta);

	% 2つの交点のうち、M3に近い方を選択し、M1とする関数
	[H1] = judgePointH1Function(crossH1, crossH2, M3); 


	% 不等式制約の定義 (c <= 0)

	% a,bの制約条件
	% 0 < b < sqrt(a^2-b^2) < a
	% 楕円の横軸長さaが円の外部より大きい、かつ、楕円の縦軸長さbが円の外部より小さい
  c(1) = -b + eps;           % b > 0
	c(2) = b - sqrt(a^2 - b^2) + eps; % b < f
	c(3) = sqrt(a^2 - b^2) - a + eps; % f < a


	% ミラー高さの制約条件
	% 0.8*d < H1M3 < d となるように制約を設定
	D = 12;  % カプセルの直径（mm）
	d = D / 2;  % カプセルの半径（mm）
	H1M3 = sqrt((H1(1) - M3x)^2 + (H1(2) - M3y)^2);
	% c(4) = 0.80 * d - H1M3;
	% c(5) = H1M3 - d;

	% 2fsinθ < d
	c(6) = 2 * sqrt(x(1)^2 - x(2)^2) * sin(x(3)) - d + eps;

	% thetaの制約条件 theta1 < theta < theta2
	[theta_range] = setThetaConstraint(a,b);
	c(7) = -x(3) + theta_range(1) + eps;
	c(8) = x(3) - theta_range(2) + eps;
	
	% 楕円外部にあるためには、ellipse_eq_M3が1より大きくなる（境界線上にないためにepsも引いている）
	ceq = H1M3 - d;
end


% 初期値
x0 = [3, 1, pi/2.2];  % 初期点
% x = [3, 1, pi/2.3];  % 初期点
[c, ceq] = constraint(x);   

disp("制約条件：C1");
disp(c(1));
disp("制約条件：C2");
disp(c(2));
disp("制約条件：C3");
disp(c(3));
% disp("制約条件：C4");
% disp(c(4));
% disp("制約条件：C5");
% disp(c(5));
disp("制約条件：C6");
disp(c(6));
disp("制約条件：C7");
disp(c(7));
disp("制約条件：C8");
disp(c(8));


% % fminconのオプション設定（最適化過程を表示する）
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

% % fminconの実行
[x_opt, fun] = fmincon(@objective, x0, [], [], [], [], [], [], @constraint, options);

% 結果の表示
disp('最適化された結果:');
disp(['a = ', num2str(x_opt(1))]);
disp(['b = ', num2str(x_opt(2))]);
disp(['theta = ', num2str(x_opt(3))]);
disp(['最小化された値: ', num2str(fun)]);


% 楕円と直線Lの交点を求める関数
function [crossW1, crossW2] = solveM1Intersections(a, b, f1x, f1y, theta)
  % シンボリック変数の定義
  syms Ex Ey;
  
  % 楕円の方程式
  ellipse_eq = (Ex^2 / a^2) + (Ey^2 / b^2) == 1;
  
  % 直線Lの方程式 (y = -tan(theta) * (x - f1x) + f1y)
  line_eq = -tan(theta) * Ex + tan(theta) * f1x + f1y;

  % 直線の方程式を楕円の方程式に代入して解く (楕円の式, 代入される変数, 代入する式), xについて解く
  % subsで2つの式を結びつけて、方程式を作成
  solutions = solve(subs(ellipse_eq, Ey, line_eq), Ex);
  % 解からx座標を取得
  crossWx = double(solutions);
  % 対応するy座標を計算
  crossWy = double(subs(line_eq, Ex, crossWx));

  % 交点の座標を表示
  crossW1 = [crossWx(1), crossWy(1)];
  crossW2 = [crossWx(2), crossWy(2)];
end

% 2つの交点のうち、M3に近い方を選択し、M1とする関数
function [W1] = judgePointFunction(crossW1, crossW2, M3)
  % crossW1,M3の距離
  W1M3 = sqrt((crossW1(1) - M3(1))^2 + (crossW1(2) - M3(2))^2);
  % crossW2,M3の距離
  W2M3 = sqrt((crossW2(1) - M3(1))^2 + (crossW2(2) - M3(2))^2);

  % crossL1とL2のうち、M3に近い方を選択し、M1とする
  if W1M3 < W2M3
      W1 = crossW1;
  else
      W1 = crossW2;
  end

end



  function [f1x, f1y, f2x, f2y, f] = solveFocus(a, b)
    % 焦点座標 f1：光源, f2：集光点
    f1x = -sqrt(abs(a^2 - b^2));
    f1y = 0;
    f2x = sqrt(abs(a^2 - b^2));
    f2y = 0;
    f = abs(f1x - f2x);
  end


function [M3x, M3y, M3] = solveM3(f1x, f2x, theta)
    % M3の座標
    F = abs(f1x - f2x)/2;
    M3x = F * cos(2 * theta);
    M3y = -F * sin(2 * theta);
    M3 = [M3x, M3y];
end




% 楕円と直線Nの交点を求める関数
function [crossH1, crossH2] = solveH1Intersections(a, b, M3x, M3y, theta)
  % シンボリック変数の定義
  syms Ex Ey;
  
  % 楕円の方程式
  ellipse_eq = (Ex^2 / a^2) + (Ey^2 / b^2) == 1;
  
  % 直線Lの方程式 (y = -tan(theta) * (x - f1x) + f1y)
  line_eq = 1 / tan(theta) * (Ex - M3x) + M3y;
 
  % 直線の方程式を楕円の方程式に代入して解く (楕円の式, 代入される変数, 代入する式), xについて解く
  % subsで2つの式を結びつけて、方程式を作成
 
  solutions = solve(subs(ellipse_eq, Ey, line_eq), Ex);
  % 解からx座標を取得
  crossHx = double(solutions);

  % 対応するy座標を計算
  crossHy = double(subs(line_eq, Ex, crossHx));

  % 交点の座標を表示
  crossH1 = [crossHx(1), crossHy(1)];
  crossH2 = [crossHx(2), crossHy(2)];
  
end

% 2つの交点のうち、M3に近い方を選択し、H1とする関数
function [H1] = judgePointH1Function(crossH1, crossH2, M3)
  % crossW1,M3の距離
  H1M3 = sqrt((crossH1(1) - M3(1))^2 + (crossH1(2) - M3(2))^2);
  % crossW2,M3の距離
  H2M3 = sqrt((crossH2(1) - M3(1))^2 + (crossH2(2) - M3(2))^2);

  % crossL1とL2のうち、M3に近い方を選択し、M1とする
  if H1M3 < H2M3
      H1 = crossH1;
  else
      H1 = crossH2;
  end

end

% thetaの制約条件を設定
function [theta_range] = setThetaConstraint(a,b)
	% theta範囲は楕円と円の交点区間で決まるため2つの交点C1,C2を求める
	% シンボリック変数の定義
	syms Ex Ey theta;

	% 楕円の方程式 (Ex^2 / a^2) + (Ey^2 / b^2) = 1
	ellipse_eq = (Ex^2 / a^2) + (Ey^2 / b^2) == 1;

	% 焦点間の距離 f (楕円の焦点距離 f = sqrt(a^2 - b^2))
	f = sqrt(a^2 - b^2);
	% 円の方程式 Ex^2 + Ey^2 = f^2
	circle_eq = Ex^2 + Ey^2 == f^2;

	% 楕円の方程式に円の方程式を代入して、xについて解く
	% Ey^2 を circle_eq から代入し、Exについて解を得る
	crossCx = solve(subs(ellipse_eq, Ey^2, f^2 - Ex^2), Ex);

	% crossCx は複数解が存在するので、それぞれの解に対して Ey を計算
	crossCy = double(sqrt(f^2 - double(crossCx).^2));

	% 2つの交点 C1 (Ex1, Ey1) と C2 (Ex2, Ey2) が得られるので、それぞれの角度を計算
	C1 = [crossCx(1), crossCy(1)];


	C1ABS = [abs(C1(1)), abs(C1(2))];
    
  XposPoint = double([C1ABS(1), -C1ABS(2)]);
  XnegPoint = double([-C1ABS(1), -C1ABS(2)]);


	% 直線Lの式から角度を求める
	f1x = -sqrt(a^2 - b^2);

	% 直線の式(pos用) y = -tan(theta) * (x - f1x)
	line_L_pos = XposPoint(2) == -tan(theta) * (XposPoint(1) - f1x);

	% 直線の式(neg用) y = -tan(theta) * (x - f1x)
	line_L_neg = XnegPoint(2) == -tan(theta) * (XnegPoint(1) - f1x);

	
	% theta について解く
	start_theta = double(solve(line_L_pos, theta));
	
	end_theta = double(solve(line_L_neg, theta));

	theta_range = [start_theta, end_theta];

end


