% MATLAB hook to call the music_doa executable and return matrix A
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

  X = csvread('output/X.csv'); # Read ULA response, one column per snapshot

  data = csvread('output/S_MUSIC.csv');
  theta = data(:, 1); % First column is angle
  S_MUSIC = data(:, 2); % Second column is spectrum

  data = csvread('output/DOA.csv');
  T = data(:, 1); % First column is true angles
  T_EST = data(:, 2); % Second column is estimated angles  

  T = sort(T)'
  T_EST = sort(T_EST)'

  figure(1);
  subplot(2, 1, 1);
%    surf(abs(X));
%    scatter(0:M-1, real(X(:,1)), 'x');
    plot(0:M-1, real(X(:,50)));
    title('ULA Snapshot');

  subplot(2, 1, 2);
    plot(theta, S_MUSIC);
    title('MUSIC Pseudo-spectrum');
    xlabel('Angle (degrees)'); ylabel('Spectrum');
    grid on;

    hold on;
%    for k = 1:length(T)
%      line([T_EST(k) T_EST(k)], ylim, 'Color', 'b', 'LineStyle', '--');
%    end

    for k = 1:length(T)
      line([T(k) T(k)], ylim, 'Color', 'k', 'LineStyle', '--');
    end

end
