% Link budget analysis and data storage
results = struct('Time', [], 'Satellite', [], 'Distance_km', [], 'LinkMargin_dB', []);

for t = times
    sc.advance(t); % Advance simulation to the next timestep
    current_time = datetime('now') + seconds(t);

    % Loop through each satellite to check line-of-sight and link margin
    % for i = 1:n_satellites
        % sat = satellites(i);
        [isLOS, distance] = accessStatus(sat_lunar1, groundStation); % Line-of-sight and distance to GS

        if isLOS
            % Calculate free-space path loss (FSPL)
            distance_km = distance / 1e3; % Convert to km
            FSPL = fspl(distance, f); % Use built-in fspl function (Antenna Toolbox)

            % Calculate received power in dBW
            Pr = Tx_power + Gt + Gr - FSPL;

            % Link margin
            noise_floor = -174 + 10*log10(f); % Simplified noise floor in dBm
            link_margin = Pr - noise_floor; % dB

            % Store results
            results.Time = [results.Time; current_time];
            results.Satellite = [results.Satellite; i];
            results.Distance_km = [results.Distance_km; distance_km];
            results.LinkMargin_dB = [results.LinkMargin_dB; link_margin];

            % Display results
            fprintf('Time: %s, Satellite: %d, Distance: %.2f km, Link Margin: %.2f dB\n', ...
                    datestr(current_time), i, distance_km, link_margin);
        end
    % end
end

%%
% % Link budget analysis and data storage
% results = struct('Time', [], 'Satellite', [], 'Distance_km', [], 'LinkMargin_dB', []);
%
% for t = times
%     sc.Advance(t); % Advance simulation to the next timestep
%     current_time = datetime('now') + seconds(t);
% 
%     % Loop through each satellite to check line-of-sight and link margin
%     for i = 1:n_satellites
%         sat = satellites(i);
%         [isLOS, distance] = accessStatus(sat, groundStation); % Line-of-sight and distance to GS
% 
%         if isLOS
%             % Calculate free-space path loss (FSPL)
%             distance_km = distance / 1e3; % Convert to km
%             FSPL = fspl(distance, f); % Use built-in fspl function (Antenna Toolbox)
% 
%             % Calculate received power in dBW
%             Pr = Tx_power + Gt + Gr - FSPL;
% 
%             % Link margin
%             noise_floor = -174 + 10*log10(f); % Simplified noise floor in dBm
%             link_margin = Pr - noise_floor; % dB
% 
%             % Store results
%             results.Time = [results.Time; current_time];
%             results.Satellite = [results.Satellite; i];
%             results.Distance_km = [results.Distance_km; distance_km];
%             results.LinkMargin_dB = [results.LinkMargin_dB; link_margin];
% 
%             % Display results
%             fprintf('Time: %s, Satellite: %d, Distance: %.2f km, Link Margin: %.2f dB\n', ...
%                     datestr(current_time), i, distance_km, link_margin);
%         end
%     end
% end

% % Plot link margin over time for each satellite
% figure;
% hold on;
% for i = 1:n_satellites
%     satellite_data = results(results.Satellite == i, :);
%     plot(satellite_data.Time, satellite_data.LinkMargin_dB, 'DisplayName', ['Satellite ' num2str(i)]);
% end
% xlabel('Time');
% ylabel('Link Margin (dB)');
% title('Link Margin Over Time for Lunar Constellation');
% legend show;
% grid on;
