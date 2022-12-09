function [tyInfo]  = plotYData(yStoreNoise, stationNum, t)    
    for i = 2:length(yStoreNoise)
        for j = 1:12
            if isempty(find(stationNum(j,i)==j)) == 1
                continue;
            else
                if isempty(yStoreNoise{i}) == 1
                    continue;
                else
                    if length(yStoreNoise{i}) == 3
                        tyInfo(:,i,find(stationNum(:,i)~=0)) = [yStoreNoise{i}; j];     
                    else
                        station_num_vec = find(stationNum(:,i) ~= 0);
                        if station_num_vec(1) ~= 1 && station_num_vec(1) ~= 12
                            tyInfo(:,i,station_num_vec(1)) = [yStoreNoise{i}(1:3); station_num_vec(1)];
                            tyInfo(:,i,station_num_vec(2)) = [yStoreNoise{i}(4:6); station_num_vec(2)];
                        else
                            tyInfo(:,i,station_num_vec(2)) = [yStoreNoise{i}(1:3); station_num_vec(2)];
                            tyInfo(:,i,station_num_vec(1)) = [yStoreNoise{i}(4:6); station_num_vec(1)];
                        end
                    end
                end
            end
        end
    end
    tyInfo(tyInfo==0) = NaN;
    plotvec1 = tyInfo(:,:,1);
    plotvec2 = tyInfo(:,:,2);
    plotvec3 = tyInfo(:,:,3);
    plotvec4 = tyInfo(:,:,4);
    plotvec5 = tyInfo(:,:,5);
    plotvec6 = tyInfo(:,:,6);
    plotvec7 = tyInfo(:,:,7);
    plotvec8 = tyInfo(:,:,8);
    plotvec9 = tyInfo(:,:,9);
    plotvec10 = tyInfo(:,:,10);
    plotvec11 = tyInfo(:,:,11);
    plotvec12 = tyInfo(:,:,12);
    
    figure
    hold on
    grid on
    plot(t,plotvec1(1,:),t,plotvec2(1,:),t,plotvec3(1,:),t,plotvec4(1,:),t,plotvec5(1,:),t,plotvec6(1,:),t,plotvec7(1,:),t,plotvec8(1,:),t,plotvec9(1,:),t,plotvec10(1,:),t,plotvec11(1,:),t,plotvec12(1,:))
    title('Noisy Range and Time')
    xlabel('Time [seconds]')
    ylabel('Range [km]')
    
    
    figure
    hold on
    grid on
    plot(t,plotvec1(2,:),t,plotvec2(2,:),t,plotvec3(2,:),t,plotvec4(2,:),t,plotvec5(2,:),t,plotvec6(2,:),t,plotvec7(2,:),t,plotvec8(2,:),t,plotvec9(2,:),t,plotvec10(2,:),t,plotvec11(2,:),t,plotvec12(2,:))
    xlabel('Time (seconds)')
    title('Noisy Range Rate and Time')
    ylabel('Range Rate [km/s]')
    
    figure
    hold on
    grid on
    plot(t,plotvec1(3,:),t,plotvec2(3,:),t,plotvec3(3,:),t,plotvec4(3,:),t,plotvec5(3,:),t,plotvec6(3,:),t,plotvec7(3,:),t,plotvec8(3,:),t,plotvec9(3,:),t,plotvec10(3,:),t,plotvec11(3,:),t,plotvec12(3,:))
    xlabel('Time [sec]')
    title('Noisy Elevation Angle and Time')
    ylabel('Elevation Angle [rad]')
    
    
end