function pic = field2pic(field)
H = (angle(field)+pi)/2/pi; % Hue
S = ones(size(field)); % Saturation
V = funcs.nmlz(abs(field)); % Value
hsv_image(:,:,1) = H;
hsv_image(:,:,2) = S;
hsv_image(:,:,3) = V;

pic = hsv2rgb(hsv_image);