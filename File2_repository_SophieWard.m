function [Uq,Vq,Hq]=unstruc_interp_slw(x1,y1,utmE,utmN,UU,VV,HH,IKLE,NT,roi_size,x_offset,y_offset,in_mask_check,hFig)

% Unstrucuted ptm, identify nearest nodes, which triangle it sits in and do
% barycentric interpolation to that point, outputs H/U/V at that point

% % % troubleshooting
% x1=xnew;  	% Point of interest
% y1=ynew; 
% utmE;     	% lat/lon in utm of the gridded parameters
% utmN;
% UU=U;     	% u/v-vels and depths of the gridded parameters
% VV=V;
% HH=H;
% IKLE=ikle; 	% Connectivity matrix 
% NT=NT;   	% Timestep of the loop / parameter
% hnew=griddata(utmE,utmN,H(:,tide),xnew,ynew) %interpolate to check if in mask

while true
    % Region of interest
    roi_x = x1 - (x_offset * roi_size / 2);
    roi_y = y1 - (y_offset * roi_size / 2);
    roi_w = x_offset * roi_size;
    roi_h = y_offset * roi_size;

	% Nodes within region of interest
    loc = find((utmE > roi_x) & (utmE < (roi_x + roi_w)) & (utmN > roi_y) & (utmN < (roi_y + roi_h)));

    % Zoom in on region of interest
%    set(hFig.Children, 'Xlim', [roi_x - 1000 roi_x + roi_w + 1000]);
%    set(hFig.Children, 'Ylim', [roi_y - 1000 roi_y + roi_h + 1000]);
%     plot(utmE(loc), utmN(loc), 'go');
            
    % Find the nearest nodes and their connected neighbours
    search_range=9;
    if numel(loc) < search_range
        roi_size=roi_size*1.5; % Increase the ROI search area until we have N nodes
        continue
    else
        % Calculate the euclidean distances for those nodes in the ROI
        roi = cat(3, utmE(loc), utmN(loc), loc);
		dist2 = sum(([x1 y1] - roi(:, 1:2)) .^ 2, 2);

		particle_in_mask = ismember(0,dist2);

        if particle_in_mask && ~in_mask_check
            % Use the current node and two nearest nodes to create
            % a triangle by the particle's coordinates.

            [~, ind] = mink(dist2, 3);             % Find the nearest nodes to the particle

            i1 = roi(ind(1),3);
            i2 = roi(ind(2),3);
            i3 = roi(ind(3),3);
            
            ax = roi(ind(1),1); ay = roi(ind(1),2);
            bx = roi(ind(2),1); by = roi(ind(2),2);
            cx = roi(ind(3),1); cy = roi(ind(3),2);

            line([ax, bx], [ay, by], 'Color', 'b');
            line([bx, cx], [by, cy], 'Color', 'b');
            line([cx, ax], [cy, ay], 'Color', 'b');
            break
                    
         else
            % Find the nearest node to the particle
            [~, ind] = mink(dist2, 1);
            [roi_nodes_ind, ~] = find((IKLE == roi(ind, 3)));

            % Nodes connected to the nearest node
            roi_nodes = IKLE(roi_nodes_ind, :);
            
            if size(roi_nodes,1)<3
                	Uq=0;
                    Vq=0;
                    Hq=0;
                    break
            else
                
                % Create a triangle around the particle's coordinates
                for i=1:length(roi_nodes)
                    ax = utmE(roi_nodes(i,1)); ay = utmN(roi_nodes(i,1));
                    bx = utmE(roi_nodes(i,2)); by = utmN(roi_nodes(i,2));
                    cx = utmE(roi_nodes(i,3)); cy = utmN(roi_nodes(i,3));

                    side_1 = (x1 - bx) * (ay - by) - (ax - bx) * (y1 - by);
                    side_2 = (x1 - cx) * (by - cy) - (bx - cx) * (y1 - cy);
                    side_3 = (x1 - ax) * (cy - ay) - (cx - ax) * (y1 - ay);

                     line([ax, bx], [ay, by], 'Color', 'b');
                     line([bx, cx], [by, cy], 'Color', 'b');
                     line([cx, ax], [cy, ay], 'Color', 'b');

                    % All signs must be positive or negative to be inside, this doesn't
                    % work if in the mask, XXX
                    if (side_1>0 && side_2>0 && side_3>0) || (side_1<0 && side_2<0 && side_3<0)
                        i1=roi_nodes(i,1);
                        i2=roi_nodes(i,2);
                        i3=roi_nodes(i,3);
                        break;
                    end    
                end
            
			end
			break
        end
    end
end

if ~exist('i1')
	Uq=0;
	Vq=0;
	Hq=0;
else
	u1=UU(i1,NT); v1=VV(i1,NT); h1=HH(i1,NT);
	u2=UU(i2,NT); v2=VV(i2,NT); h2=HH(i2,NT);
	u3=UU(i3,NT); v3=VV(i3,NT); h3=HH(i3,NT);

	m1=[utmE(i1) utmN(i1)];
	m2=[utmE(i2) utmN(i2)];
	m3=[utmE(i3) utmN(i3)];

	atri=abs((m1(1)*m2(2) - m1(1)*m3(2) + m2(1)*m3(2) - m2(1)*m1(2) + m3(1)*m1(2) - m3(1)*m2(2))/2); %area(m1m2m3)
	a1=abs((m2(1)*y1 - m2(1)*m3(2) + x1*m3(2) - x1*m2(2) + m3(1)*m2(2) - m3(1)*y1)/2);% ./ atri %area(m2pm3)
	a2=abs((m1(1)*y1 - m1(1)*m3(2) + x1*m3(2) - x1*m1(2) + m3(1)*m1(2) - m3(1)*y1)/2);% ./ atri %area(m1pm3)
	a3=abs((m1(1)*y1 - m1(1)*m2(2) + x1*m2(2) - x1*m1(2) + m2(1)*m1(2) - m2(1)*y1)/2);% ./ atri %area(m1pn2)

	% Calculate coefficients
	alpha1=a1/atri; alpha2=a2/atri; alpha3=a3/atri;

	% Calculate velocity at point of interest,
	Uq=alpha1*u1 + alpha2*u2 + alpha3*u3;
	Vq=alpha1*v1 + alpha2*v2 + alpha3*v3;
	Hq=alpha1*h1 + alpha2*h2 + alpha3*h3;
    
    clear i1 i2 i3
end
end
