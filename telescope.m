%% Correction to telescope observation

%% Input the pointed ra (deg), dec (deg), the velocity vector of the
%% satellite at the expected time, and the angle between the entrance
%% and exit from the field of view.

%% Returns the right ascension and declination in degrees of the point
%% of closest approach to the center of the field of view.

%% Presumes the spotting scope on the 14" schmidt-cassegrain at the UM
%% observatory with a 123" field of view; if not, correct fieldofview
%% accordingly.

%% L. Healy, ENAE 441, Oct. 4, 2008.

%% format long g
%% telescope(310, -28, [6.691082315278334, 3.289744906327188, 1.823770927897560], 74)
%% ans =
%%         -50.2946717586475         -27.2241575246794

%%function [ra, dec] = telescope(ra, dec, velocity, anglebetween)
function vec = telescope(ra, dec, velocity, anglebetween)
  fieldofview = 123/60*pi/180*1/2;
  rad = ra*pi/180;
  decd = dec*pi/180;

  obsdir = [cos(rad)*cos(decd), sin(rad)*cos(decd), sin(decd)];
  
  vhat = velocity - obsdir * dot(obsdir,velocity);
  vhatp = vhat/norm(vhat);
  rotax = vhatp;
  ang = -fieldofview*cos(anglebetween*pi/360);

  %% Do an active 3D rotation about rotax by an angle ang
  cosang = cos(-ang);
  rot...
  = cosang*eye(3,3)...
  + (1-cosang)*[rotax(1)*rotax(1), rotax(2)*rotax(1), rotax(3)*rotax(1);
		rotax(1)*rotax(2), rotax(2)*rotax(2), rotax(3)*rotax(2);
		rotax(1)*rotax(3), rotax(2)*rotax(3), rotax(3)*rotax(3)]...
  + sin(-ang)*[0, rotax(3), -rotax(2) ;
	      -rotax(3), 0, rotax(1);
	      rotax(2), -rotax(1), 0];

  new = rot*obsdir';
  
  180/pi*[atan2(new(2),new(1)), atan(new(3)/norm(new(1:2)))]
