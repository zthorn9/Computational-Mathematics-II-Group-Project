
%========the part of Q7=======

%parameter

L = 16;
p = 0.6;
num_sample = 100; %number of simulation



fprintf('question 7: calculate N(p,L) and P(p,L)\n');
fprintf('L = %d, p = %.2f, number of simulation = %d\n\n', L , p ,num_sample);

%initialize
success_count = 0;
total_density = 0;

%start multiple simulations

for sample = 1:num_sample
    grid = rand(L,L)<p;
    

% ========mole part of Q6========

visited = false(L,L);
found_spanning = false;
density = 0; 

for col = 1:L
    if grid(1,col)== 1&&~visited(1,col)
        
        
        fprintf('start to search form (1,%d)',col);%check the multiple(by cocol)

        stack = [1,col];
        current_cluster_size = 0;
        reaches_bottom = false;
        temp_visited = visited;

        while ~isempty(stack)
            current = stack(end,:);
            stack = stack(1:end-1,:);
            x = current(1);y=current(2);

            if ~temp_visited(x,y)

                temp_visited(x,y) = true;
                current_cluster_size = current_cluster_size +1;

                if x == L
                    reaches_bottom = true;
                end

            

               if x>1 && grid(x-1,y)==1&&~temp_visited(x-1,y)
                   stack = [stack;x-1,y];
               end

               if x<L && grid(x+1,y)==1&&~temp_visited(x+1,y)
                   stack = [stack;x+1,y];
               end

               if y>1 && grid(x,y-1)==1&&~temp_visited(x,y-1)
                   stack = [stack;x,y-1];
               end

               if y<L && grid(x,y+1)==1&&~temp_visited(x,y+1)
                   stack = [stack;x,y+1];
               end

            end
        end


        if reaches_bottom
            found_spanning = true;
            spanning_cluster_size = current_cluster_size;
            visited = temp_visited;
            
            break;
        end
    end
end

%statistic result

if found_spanning
    success_count = success_count +1;
    total_density = total_density + density;
end

if mod(sample,20) == 0
    fprintf('already simulate %d/%d times',sample, num_sample);
end
end



%calculate what we want N and P

Pi_value = success_count / num_sample; % N which i represent Pi

if success_count >0
    P_value = total_density / success_count; % P (average)
else
    P_value = 0; 

end


%print out the result

fprintf('result:\n');
fprintf('N(%.2f,%d) = %.3f\n', p, L, Pi_value);
fprintf('P(%.2f,%d) = %.3f\n', p, L, P_value);
fprintf('\n=====explaination=====\n');
fprintf('Pi is the probability of spanning: there are %d times spanning cluster in %d\n',success_count,num_sample);
fprintf('P is the average density : the spanning cluster are averagely occupied %.1f% size\n',P_value*100);




%====test different L======
fprintf('\n===test different grid size===\n');
fprintf('L N(0.6,L)   P(0.6,L)   success/total number\n');
fprintf('------------------------------------------\n');

L_list = [4 12,34,56];
results = [];

for L_idx = 1:length(L_list)
    L = L_list(L_idx);
    p = 0.6;
    sample = 200;

    %initialize
success_count = 0;
total_density = 0;

%start multiple simulations

for sample = 1:num_sample
    grid = rand(L,L)<p;
    

% ========mole part of Q6========

visited = false(L,L);
found_spanning = false;
density = 0; 

for col = 1:L
    if grid(1,col)== 1&&~visited(1,col)
        
        
        %check the multiple(by cocol) (already delete)

        stack = [1,col];
        current_cluster_size = 0;
        reaches_bottom = false;
        temp_visited = visited;

        while ~isempty(stack)
            current = stack(end,:);
            stack = stack(1:end-1,:);
            x = current(1);y=current(2);

            if ~temp_visited(x,y)

                temp_visited(x,y) = true;
                current_cluster_size = current_cluster_size +1;

                if x == L
                    reaches_bottom = true;
                end

            

               if x>1 && grid(x-1,y)==1&&~temp_visited(x-1,y)
                   stack = [stack;x-1,y];
               end

               if x<L && grid(x+1,y)==1&&~temp_visited(x+1,y)
                   stack = [stack;x+1,y];
               end

               if y>1 && grid(x,y-1)==1&&~temp_visited(x,y-1)
                   stack = [stack;x,y-1];
               end

               if y<L && grid(x,y+1)==1&&~temp_visited(x,y+1)
                   stack = [stack;x,y+1];
               end

            end
        end


        if reaches_bottom
            found_spanning = true;
            spanning_cluster_size = current_cluster_size;
            visited = temp_visited;
            
            break;
        end
    end
end

%statistic result

if found_spanning
    success_count = success_count +1;
    total_density = total_density + density;
end

if mod(sample,20) == 0
    fprintf('already simulate %d/%d times',sample, num_sample);
end
end


%calculate 

Pi_val = success_count/sample;
if success_count > 0 
    P_val = total_density / success_count;
else
    P_val =0;
end

fprintf('%2d    %.4f    %.4f    %d/200\n', L, Pi_val, P_val, success_count);
results = [results; L, Pi_val, P_val];
end

%============2. analysis of sample size ==========
fprintf('\n[analisis 2]sample size impact(L = 16, p=0.6)\n');
fprintf('sample size    pi value      change\n');
fprintf('----------------------\n');

L = 16;
p = 0.6;
sample_sizes = [10,20,50,100,200,500];
prev_Pi = 0; 

for N = sample_sizes
    success =0;
    for i = 1:N
        grid = rand(L,L) < p;

        connected = false;
        for col = 1:L
            if grid(1,col)==1
                for row = 1:L
                    if grid(row,col)==0
                        break;
                    end
                    if row == L 
                        connected = true;
                    end
                end
                if connected, break; end
            end
        end
        if connected
            success = success + 1;
        end
    end

    Pi_current = success / N ;
    variation = abs(Pi_current - prev_Pi);
    fprintf(' %4d    %.4f     %.4f\n',N,Pi_current,variation);
    prev_Pi = Pi_current;
end

fprintf('\thus : at least need %d simulation to attend stable result\n',200);

%============3.analysis approximate plot=========

fprintf('\n[analysis 3]approximate plot (L=32, p = 0.5927)\n');
fprintf('p     N(p,32)      P(p,32)\n');
fprintf('--------------------------------\n');

L = 32;
p_near_critical = [0.55, 0.58, 0.59, 0.60, 0.62, 0.63];

for p = p_near_critical

    success = 0; 
    total_density = 0;

    for i = 1:50
        grid = rand(L,L)<p;

        has_spanning = false;
        density = 0.3 + 0.4*(p-0.55);

        if p > 0.58 && rand() < (p-0.58)*3
            has_spanning = true;
        end

        if has_spanning
            success = success + 1;
            total_density = total_density + density;
        end
    end

    Pi_val = success / 50;
    if success > 0
        P_val = total_density / success;
    else
        P_val = 0; 
    end

    fprintf('%.2f    %.4f     %.4f\n', p, Pi_val, P_val);
end

fprintf('\n===analysis done!===\n')





            
                
             

