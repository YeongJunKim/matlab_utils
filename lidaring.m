%% KOREA UNIVERCITY
%% author : colson@korea.ac.kr(dud3722000@naver.com)
%% For D-FIR (lidar clustering)

classdef lidaring < handle
    properties
        %% variable area
        angle_interval  = 0;
        time_interval   = 0;
        robot_num       = 0;
        lidar_range     = [];
        
        %% robot spec
        % linear and angular;
        maximum_velocity = [];
        %% get data
        scan_data = [];
        interested_scan = [];
        detected_num = [];
        %% result data
        history_size = 0;
        result_Cartesian = [];
        result_Cluster = [];
        result = [];
        %% robot position
        result_data = [];
        result_mean = [];
        result_now = [];
        
        step = 0;
        
        %% class info
        namespace = "";
        fig;
        is_init = "no";
    end
    methods
        %% init
        function obj = lidaring_init(obj, namespace, angle_interval, time_interval, robot_num, range)
            obj.namespace = strcat(namespace,"/scan");
            obj.fig             = figure('Name', obj.namespace);
            obj.angle_interval  = angle_interval;
            obj.time_interval   = time_interval;
            obj.robot_num       = robot_num;
            obj.lidar_range     = range;
            
            factor              = 360/angle_interval;
            obj.maximum_velocity = zeros(2, obj.robot_num);
            obj.scan_data       = zeros(1, factor);
            
            obj.is_init = "ok";
        end
        function r = lidaring_position_init(obj, init, size)
            obj.history_size = size;
            obj.step = 1;
            obj.result_data = zeros(2, obj.robot_num, size);
            obj.result_data(:,:,1) = init;
            r = "ok";
        end
        %% run
        function r = lidaring_run(obj, time_interval, scan_data, tracking_position, theta)
            if(obj.step < obj.history_size)
                obj.step = obj.step+1;
                obj.scan_data = scan_data;
                scan = lidarScan(scan_data, deg2rad(1:360));
                obj.interested_scan = removeInvalidData(scan, 'RangeLimits', obj.lidar_range);
                %             disp(interested_scan.Count)
                obj.detected_num = 0;
                if(obj.interested_scan.Count == 0)
                    disp("[lidaring.m] No detection - alternate before data")
                    
                    obj.result_data(:,:,obj.step) = obj.result_data(:,:,obj.step - 1);
                    
                elseif(obj.interested_scan.Count == 1)
                    disp("[lidaring.m] One object detected ")
                    obj.detected_num = 1;
                    obj.result_Cartesian = obj.interested_scan.Cartesian;
                    obj.result_Cluster = [1];
                    obj.result = [obj.result_Cartesian obj.result_Cluster];
                else
                    disp("[lidaring.m] More than one objects detected ")
                    Z = linkage(obj.interested_scan.Cartesian, 'single');
                    C = cluster(Z, 'cutoff', 0.15, 'criterion','distance');
                    
                    obj.detected_num = max(C);
                    obj.result_Cartesian = obj.interested_scan.Cartesian;
                    obj.result_Cluster = C;
                    obj.result = [obj.result_Cartesian obj.result_Cluster];
                end
                obj.result_mean = zeros(2, obj.detected_num);
                
                
                
                for i = 1:obj.detected_num
                    k = find(obj.result(:,3)==i);
                    obj.result_mean(:, i) = sum(obj.result_Cartesian(k,:))./size(k,1);
                end
                
                
                distance_factor = zeros(1, obj.detected_num);
                for i = 1:obj.robot_num
                    for j = 1:obj.detected_num
                        x_diff = obj.result_mean(1,j) - obj.result_data(1, i,obj.step - 1);
                        y_diff = obj.result_mean(2,j) - obj.result_data(2, i,obj.step - 1);
                        distance_factor(1,j) = x_diff^2 + y_diff^2;
                    end
%                     disp(obj.namespace);
                    a = zeros(1, size(distance_factor,2));
                    a(:) = distance_factor(1,:);
                    [M, I] = min(distance_factor(:));
                    [row, col] = ind2sub(size(distance_factor), I);
%                     who_is = zeros(size(distance_factor,2),1);
%                     who_is(col) = 1;
%                     obj.result_mean;
%                     a;
%                     who_is;
%                     point_with_distance_factor = [obj.result_mean' a' who_is];
                    obj.result_data(:,i,obj.step) = obj.result_mean(:,col);
                    obj.result_now = obj.result_mean(:,col);
                    
                    if norm(obj.result_data(:,i,obj.step)) < 0.15
                    obj.result_data(:,i,obj.step) = obj.result_data(:,i,obj.step-1);
                    obj.result_now = obj.result_data(:,i,obj.step);
                    end
                        
                    dist = norm(obj.result_data(:,i,obj.step) - obj.result_data(:,i,obj.step -1));
                    if dist > time_interval * 0.22 * 4
                       disp("[lidaring.m] Scan error");
                       obj.result_data(:,i,obj.step) = obj.result_data(:,i,obj.step -1);
                    end
                end
                
                
                
                
                %% Based on known robot dynamics expected position is?
                
                if(obj.namespace == "/tb3a/sca1n")
                    lidaring_draw_now_prediction(obj, 1, 0.22 * obj.time_interval);
                elseif(obj.namespace == "/tb3b/sc1an")
                    lidaring_draw_now_prediction(obj, 1, 0.22 * obj.time_interval);
                elseif(obj.namespace == "/tb3c/s1can")
                    lidaring_draw_now_prediction(obj, 1, 0.22 * obj.time_interval);
                end
            else
                disp("[lidaring.m] Lidar over history");
            end
            
            
            
            
            r = "ok";
        end
        function r = lidaring_draw_now_prediction(obj, robot_num, radius)
            figure(obj.fig);
            clf;
            plot(obj.interested_scan);
            hold on;
            %             obj.result_data(:,robot_num, obj.step)
            plot(obj.result_mean(1,:), obj.result_mean(2,:), "*");
            hold on;
            temp = zeros(1,obj.robot_num);
            temp(:) = radius;
            viscircles(obj.result_data(:, :, obj.step)', temp);
            
            xlim([-1 1]);
            ylim([-1 1]);
            grid on;
        end
        function r = lidaring_draw_circle(obj, radius)
            %             figure(obj.fig);
            clf;
            plot(obj.interested_scan);
            hold on;
            plot(obj.result_mean(1,:), obj.result_mean(2,:), "*");
            hold on;
            rad = zeros(size(obj.result_mean,2), 1);
            rad(:) = radius;
            viscircles(obj.result_mean',rad);
            xlim([-1.5 1.5]);
            ylim([-1.5 1.5]);
            grid on;
        end
        function r = lidaring_mean_point(obj)
            figure(obj.fig);
            clf;
            plot(obj.interested_scan);
            hold on;
            plot(obj.result_mean(1,:), obj.result_mean(2,:), "*");
            xlim([-5 5]);
            ylim([-5 5]);
            grid on;
            drawnow;
        end
        function r = lidaring_graph(obj)
            figure(obj.fig);
            clf;
            plot(obj.interested_scan);
            xlim([-5 5]);
            ylim([-5 5]);
            grid on;
            drawnow;
        end
        function r = lidaring_scatter(obj)
            figure(obj.fig);
            clf;
            gscatter(obj.result_Cartesian(:,1), obj.result_Cartesian(:,2), obj.result_Cluster)
            xlim([-5 5]);
            ylim([-5 5]);
            grid on;
            drawnow;
        end
    end
end

%% general functions

