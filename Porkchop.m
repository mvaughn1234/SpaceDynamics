clear all; close all; clc; format long; format compact;
depart_start = ymdhms2jd(2020,4,1,12,0,0);
arrive_start = ymdhms2jd(2020,9,28,12,0,0);
u_sun = 132712440018;

depart_current = depart_start;
days_to_seconds = 24*60*60;
duration = [45:1:500];
depart_range = [0:1:200];
dv_tot_min = 1000;
for j = depart_range
    depart_current = depart_current+1;
    for i = duration
        arrive_current = depart_current+i;
        if(arrive_current >= arrive_start)
            [r_e,v_e] = findEarth(depart_current);
            [r_m,v_m] = findMars(arrive_current);
            v_en = norm(v_e); v_mn = norm(v_m);
            [tr_e_s,tr_m_s] = Lamberts_Solver(r_e,r_m,i*days_to_seconds,u_sun,1);
            tr_e_sn = norm(tr_e_s); tr_m_sn = norm(tr_m_s);
            dv_e_short = sqrt(v_en^2+tr_e_sn^2-2*tr_e_sn*v_en*(dot(v_e,tr_e_s)/(v_en*tr_e_sn)));
            dv_m_short = sqrt(v_mn^2+tr_m_sn^2-2*tr_m_sn*v_mn*(dot(v_m,tr_m_s)/(v_mn*tr_m_sn)));
            dv_tot_short = dv_e_short+dv_m_short;
            [tr_e_l,tr_m_l] = Lamberts_Solver(r_e,r_m,i*days_to_seconds,u_sun,-1);
            tr_e_ln = norm(tr_e_l); tr_m_ln = norm(tr_m_l);
            dv_e_long = sqrt(v_en^2+tr_e_ln^2-2*tr_e_ln*v_en*(dot(v_e,tr_e_l)/(v_en*tr_e_ln)));
            dv_m_long = sqrt(v_mn^2+tr_m_ln^2-2*tr_m_ln*v_mn*(dot(v_m,tr_m_l)/(v_mn*tr_m_ln)));
            dv_tot_long = dv_e_long+dv_m_long;
            if(dv_tot_short <= dv_tot_long)
                c3 = dv_e_short^2;
                vinf = dv_m_short;
                if(dv_tot_short < dv_tot_min)
                    dv_tot_min = dv_tot_short;
                    min_dv_duration = i;
                    min_depart_date = depart_current;
                    min_arrive_date = arrive_current;
                    direction = 'short';
                end
            else
                c3 = dv_e_long^2;
                vinf = dv_m_long;
                if(dv_tot_long < dv_tot_min)
                    dv_tot_min = dv_tot_long;
                    min_dv_duration = i;
                    min_depart_date = depart_current;
                    min_arrive_date = arrive_current;
                    direction = 'long';
                end
            end
            depart_index = depart_current-depart_start+1;
            arrive_index = arrive_current-arrive_start+1;
            c3_A(depart_index,arrive_index) = c3;
            vinf_A(depart_index,arrive_index) = vinf;
            tof_A(depart_index,arrive_index) = i;
        end
    end
end
figure; hold on;
contour(tof_A',[50:50:500],'black','ShowText','on');
contour(c3_A',[0:2:30,30:10:40,40:20:100],'r','ShowText','on');
contour(vinf_A',[0:1:10],'b','ShowText','on');
legend('TOF [days]','C3 [km^2/s^2]','Vinf [km/s]','Location','NorthWest');
xlabel('Days Past April 1, 2020');
ylabel('Days Past Sept 28, 2020');
hold off;