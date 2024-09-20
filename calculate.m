% カプセル
D = 12;  % カプセルの直径（mm）
d = D / 2;  % カプセルの半径（mm）


% 楕円のパラメータ
a = 3.9148;
b = 0.00036153;
theta = 0.8731;

% 楕円の焦点の計算
	f1x = -sqrt(a^2 - b^2);
	f1y = 0;
	f2x = sqrt(a^2 - b^2);
	f2y = 0;
	f1 = [f1x, f1y];
	f2 = [f2x, f2y];
	F = sqrt((f1x - f2x)^2);
	f = F / 2;

% M3の座標の計算
M3x = f * cos(2 * theta);
M3y = -f * sin(2 * theta);
M3 = [M3x, M3y];

% 点f1と点M3を通る直線の方程式 L
x = linspace(-8, 8, 100);
line_L = -tan(theta) * (x - f1x) + f1y;

% 直線Lに平行で垂直距離d離れた直線の方程式 S
line_S =  -tan(theta) * (x - f1x) + f1y + d / cos(theta);

% 点f2と点M3を通る直線の方程式 N
line_N = 1 / tan(theta) * (x - M3x) + M3y;




function [crossLX, crossLY] = solveM1Intersections(a, b, f1x, f1y, theta)
    % シンボリック変数の定義
    syms x y;
    
    % 楕円の方程式
    ellipse_eq = (x^2 / a^2) + (y^2 / b^2) == 1;
    
    % 直線の方程式 (y = -tan(theta) * (x - f1x) + f1y)
    line_eq = -tan(theta) * x + tan(theta) * f1x + f1y;

    % 直線の方程式を楕円の方程式に代入して解く (楕円の式, 代入される変数, 代入する式), xについて解く
    % subsで2つの式を結びつけて、方程式を作成
    solutions = solve(subs(ellipse_eq, y, line_eq), x);
    
    % 解からx座標を取得
    crossLX = double(solutions);
    
    % 対応するy座標を計算
    crossLY = double(subs(line_eq, x, crossLX));
end


[crossLX, crossLY] = solveM1Intersections(a, b, f1x, f1y, theta);

% 交点の座標を表示
crossL1 = [crossLX(1), crossLY(1)];
crossL2 = [crossLX(2), crossLY(2)];

% crossL1とL2のうち、M3に近い方を選択し、M1とする
if sqrt((crossL1(1) - M3(1))^2 + (crossL1(2) - M3(2))^2) < sqrt((crossL2(1) - M3(1))^2 + (crossL2(2) - M3(2))^2)
		M1x = crossL1(1);
		M1y = crossL1(2);
		M1 = [M1x, M1y];
else
		M1x = crossL2(1);
		M1y = crossL2(2);
		M1 = [M1x, M1y];
end

function [crossNX, crossNY] = solveM2Intersections(a, b, M3x, M3y, theta)
	% シンボリック変数の定義
	syms x y;
	
	% 楕円の方程式
	ellipse_eq = (x^2 / a^2) + (y^2 / b^2) == 1;
	
	% 直線の方程式 
	line_eq = 1 / tan(theta) * (x - M3x) + M3y;

	% 直線の方程式を楕円の方程式に代入して解く (楕円の式, 代入される変数, 代入する式), xについて解く
	% subsで2つの式を結びつけて、方程式を作成
	solutions = solve(subs(ellipse_eq, y, line_eq), x);
	
	% 解からx座標を取得
	crossNX = double(solutions);
	
	% 対応するy座標を計算
	crossNY = double(subs(line_eq, x, crossNX));
end


[crossNX, crossNY] = solveM2Intersections(a, b, M3x, M3y, theta);

% 交点の座標を表示
crossN1 = [crossNX(1), crossNY(1)];
crossN2 = [crossNX(2), crossNY(2)];

% crossN1とN2のうち、M3に近い方を選択し、M1とする
if sqrt((crossN1(1) - M3(1))^2 + (crossN1(2) - M3(2))^2) < sqrt((crossN2(1) - M3(1))^2 + (crossN2(2) - M3(2))^2)
		M2x = crossN1(1);
		M2y = crossN1(2);
		M2 = [M2x, M2y];
else
		M2x = crossN2(1);
		M2y = crossN2(2);
		M2 = [M2x, M2y];
end


% 直線Lの描画
plot(x, line_L, 'r--', 'DisplayName', 'Line L');
hold on;

% 直線Sの描画
plot(x, line_S, 'black', 'DisplayName', 'Line S');

% 直線Nの描画
plot(x, line_N, 'b--', 'DisplayName', 'Line N');

% 楕円の描画
t = linspace(0, 2*pi, 100);
ellipse_x = a * cos(t);
ellipse_y = b * sin(t);
plot(ellipse_x, ellipse_y, 'k-', 'DisplayName', 'Ellipse');

% 半径fの円の描画
Elitheta = linspace(0, 2 * pi, 100);
circle_x = f * cos(Elitheta);
circle_y = f * sin(Elitheta);
plot(circle_x, circle_y, 'g-', 'DisplayName', 'Circle f');






% 点 M1, M2, M3, 焦点のプロット
plot(M1(1), M1(2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Point M1');
plot(M2(1), M2(2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Point M2');
plot(M3(1), M3(2), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Point M3');
plot(f1(1), f1(2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Focus F1');
plot(f2(1), f2(2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Focus F2');

% テキスト表示
text(M1(1), M1(2), ' M1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(M2(1), M2(2), ' M2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(M3(1), M3(2), ' M3', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(f1(1), f1(2), ' f1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(f2(1), f2(2), ' f2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');



% 直線Nの式をグラフの任意の位置に表示
text(1, -6, ...
    sprintf('L：　y = -tan(\\theta) * (x - %.2f) + %.2f', f1x, f1y), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'red');

text(1, -7, ...
    sprintf('N：　y = 1 / tan(\\theta) * (x - %.2f) + %.2f', M3x, M3y), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'blue');

% 楕円の式を表示
text(1, -5, ...
    sprintf('Ellipse：　(x^2 / %.2f^2) + (y^2 / %.2f^2) = 1', a, b), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12, 'Color', 'black');




% 半径とθの範囲を設定
radius = 0.6;
theta_start = 0;  % 開始角度
theta_end = -theta;    % 終了角度

% θに基づいて円弧上の点を計算
theta = linspace(theta_start, theta_end, 100);  % θの範囲で点を生成
x_arc = f1x + radius * cos(theta);  % x座標
y_arc = f1y + radius * sin(theta);  % y座標

% 円弧の描画
plot(x_arc, y_arc, 'r', 'LineWidth', 2, 'DisplayName', 'Arc of Radius 2');









% グラフの設定
axis equal;
legend show;
xlabel('X');
ylabel('Y');
title('Ellipse and Related Points');
grid on;
hold off;




% 軸表示設定
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% 軸の設定
xlim([-6 6]);  % x 軸の範囲設定
ylim([-6 6]);  % y 軸の範囲設定
set(gca, 'XTick', -10:2:10); % x 軸の目盛りを1ずつ設定
set(gca, 'YTick', -10:2:10); % y 軸の目盛りを1ずつ設定
