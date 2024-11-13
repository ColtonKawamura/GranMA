function fftGaus(time_vector, particle_position_range, number_of_particles_to_view, x0, x_all, y_all)
% particle_position_range = [index closes to wall, index further from wall] than you want to look at

for target_distance = linspace(particle_position_range(1), particle_position_range(end), number_of_particles_to_view)
    [~, idx_particle] = min(abs(x0 - target_distance));
    figure;
    plot(time_vector(1:length(x_all(idx_particle, :))), x_all(idx_particle, :));% length of time vectory is weird because the oscillation stops before all the time finsihes
    grid on 

   % Pause the program and wait for user input for target_time
    disp('Pause... Please enter a target time value:');
    target_time = input('Enter target_time: ');  % Wait for user inpu 

    % split the time vector
    target_time = 97;
    [~, idx_time_split] = min(abs(time_vector - target_time));
    time_vector_before = time_vector(1:idx_time_split);
    time_vector_after = time_vector(idx_time_split+1:idx_time_split + length(x_all(idx_particle, idx_time_split+1:end)));
    particle_x_before = x_all(idx_particle, 1:idx_time_split);
    particle_x_after =  x_all(idx_particle, idx_time_split+1:end);
    particle_y_before = y_all(idx_particle, 1:idx_time_split);
    particle_y_after =  y_all(idx_particle, idx_time_split+1:end);

    % visually verify
    figure
    plot(time_vector_before, particle_x_before); 
    title("particle x before")
    % figure
    % plot(time_vector_before, particle_y_before);
    % title("particle y before") 
    figure;
    plot(time_vector_after, particle_x_after);
    title("particle x after")
    % figure;
    % plot(time_vector_after, particle_y_after);
    % title("particle y after")
    % fft time
    % fft parameters from the big function
    average_dt = mean(diff(time_vector_before));  % Average time step
    sampling_freq = 1 / average_dt;  % Sampling frequency
    nyquist_freq = sampling_freq / 2;  % Nyquist frequency
    freq_vector = linspace(0, 1, fix(length(time_vector)/2)+1) * nyquist_freq; % creates vector of posible frequecies from 0 to the nyquist frequency (upper limit)    grid on


    %%%%%

    % FFT for particle_x_before
    position_nn = particle_x_before - mean(particle_x_before);  % Center the data
    data = position_nn
    fps = 1 / mean(diff(time_vector_before));
    plotfft(data, fps, "b")
    title("fft  x before")
    % % FFT for particle_y_before
    % position_nn = particle_y_before - mean(particle_y_before);  % Center the data
    % data = position_nn
    % fps = 1 / mean(diff(time_vector));
    % plotfft(data, fps)
    % title("fft  y before")
    % FFT for particle_x_after
    position_nn = particle_x_after - mean(particle_x_after);  % Center the data
    data = position_nn
    fps = 1 / mean(diff(time_vector));
    plotfft(data, fps, "r")
    title("fft  x after")
    % % FFT for particle_y_after
    % position_nn = particle_y_after - mean(particle_y_after);  % Center the data
    % data = position_nn
    % fps = 1 / mean(diff(time_vector));
    % plotfft(data, fps)
    % title("fft  y after")
end

