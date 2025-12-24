% MATLAB hook to call the ./music executable and return matrix X
%
% Inputs:
%   M - number of sensors
%   N - number of snapshots
%   K - number of sources
%   T - Source angles
%   P - Source powers
%
function [X, S_MUSIC, T_EST] = music(M, N, K, T = [], P = ones(1, K), Pn = 0.05)
  % Call ./music to generate .csv outputs
  command = sprintf('./music -M %d -N %d -K %d -s %.3f', M, N, K, Pn);

  for k = 1:length(T)
    command = sprintf('%s S%d=%.2f@%.2f', command, k, T(k), P(k));
  end

  printf("%s\n", command);
  system(command);

  X = csvread('output/X.csv'); # ULA response, one column per snapshot

  data = csvread('output/S_MUSIC.csv');
  theta = data(:, 1);   % First column is angle
  S_MUSIC = data(:, 2); % Second column is spectrum

  data = csvread('output/DOA.csv');
  T = data(:, 1);     % First column is true angles
  T_EST = data(:, 2); % Second column is estimated angles  

  T = sort(T)'
  T_EST = sort(T_EST)'

  S_MUSIC = S_MUSIC / max(S_MUSIC);

  figure(1, 'Name', 'MUSIC Spectrum');

  subplot(1, 4, 1:2);
    plot(theta, S_MUSIC, 'k');
    xlabel('Angle (degrees)'); ylabel('Spectrum');
    grid on;

    hold on;
    for k = 1:length(T)
      line([T(k) T(k)], ylim, 'Color', 'b', 'LineStyle', '--');
    end

  subplot(1, 4, 3:4);
    h = polar(theta*pi/180, S_MUSIC);
    set(h, 'LineStyle', '-.');
    set(h, 'Color', 'k');
end
