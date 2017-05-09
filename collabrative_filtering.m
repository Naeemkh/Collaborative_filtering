% Machine Learning - Comp 8745 
% Project1 - Collabrative filtering based on Breese et al (1998)
% Prof. Deepak Venugopal
% Written by: Naeem Khoshnevis
% Last update: April 9, 2017

clc
clear all 
close all 

%% Load data
input_data_name='ratings50';
input_data = sprintf('%s%s',input_data_name,'.txt');
movie_rate = csvread(input_data);
% format movie,user,rating

%% Generate a user-item matrix
try 
    sname = sprintf('%s%s%s','recommend_',input_data_name,'.mat');
    load(sname)

    users=recommend.users;
    movies=recommend.movies;
    size_m=recommend.size_m;
    size_u=recommend.size_u;
    size_d=recommend.size_d;
    movie_rate_mat=recommend.movie_rate_mat;
    rated_mat=recommend.rated_mat;
    uv_bar=recommend.uv_bar;
    A1=recommend.A1;
    A2=recommend.A2;
    user_weight_mat=recommend.user_weight_mat;  
    
    
catch    
% generate a unique value of users and movies
users  = sort(unique(movie_rate(:,2)));
movies = sort(unique(movie_rate(:,1)));

% size
size_m = size(movies,1);
size_u = size(users,1);
size_d = size(movie_rate,1);

movie_rate_mat=zeros(size_u,size_m);


for i=1:size_d
    
   movie_no = find((movies==movie_rate(i,1))==1);
   movie_id = movies(movie_no,1);
   
   user_no = find((users==movie_rate(i,2))==1);
   user_id = users(user_no,1);
   
   movie_rate_mat(user_no,movie_no)=movie_rate(i,3);
    
end


% Generate rated items matrix
rated_mat = movie_rate_mat>0;

% Generate user mean vector (uv_bar(i))
uv_bar = sum(movie_rate_mat,2)./sum(rated_mat,2);

% Generate (A1 = v_aj - v_bar)
A1 = (movie_rate_mat - repmat(uv_bar,1,size_m)).*rated_mat; 

% Generate (A2 = A1.^2;)
A2 = A1.^2;

% Generate user weight matrix
user_weight_mat = zeros(size_u,size_u);



tic
for i=1:size_u
    for j=(i+1):size_u
    
        % numinator
        AA = sum(A1(i,:).*A1(j,:).*rated_mat(i,:).*rated_mat(j,:));
        % This line added to make sure that we don't store not practically
        % usefule data.
%         min_ncm = 3; % min number of common movies
%         if sum(rated_mat(i,:).*rated_mat(j,:)) < min_ncm
%            AA = 0; 
%         end
        
        % denuminator
        BB = sqrt(sum(A2(i,:).*rated_mat(i,:).*rated_mat(j,:)).*sum(A2(j,:).*rated_mat(i,:).*rated_mat(j,:)));
        
        % Weight_a,i
        if (BB ~= 0)
        user_weight_mat(i,j)=AA./BB;
        user_weight_mat(j,i)=AA./BB;
        else
        user_weight_mat(i,j)=0;
        user_weight_mat(j,i)=0;
        end
        
    end
end
user_weight_time=toc;
% store the data in structure
recommend.users=users;
recommend.movies=movies;
recommend.size_m=size_m;
recommend.size_u=size_u;
recommend.size_d=size_d;
recommend.movie_rate_mat=movie_rate_mat;
recommend.rated_mat=rated_mat;
recommend.uv_bar=uv_bar;
recommend.A1=A1;
recommend.A2=A2;
recommend.user_weight_mat=user_weight_mat;
recommend.user_weight_time=user_weight_time;

% save the structure
sname = sprintf('%s%s%s','recommend_',input_data_name,'.mat');
save(sname,'recommend','-v7.3');

end
% toc
stop_flag=0;

disp('Please enter user id and number of recommended movie with comma separation.')
disp('Enter 0,0 to exit')
while stop_flag==0 


input_value = input('Please enter (UserID,K): ','s'); 
out=regexp(input_value,',','split');
try
i=str2double(out{1});
k=str2double(out{2});
catch
disp('The format is not correct. Please try again: UserID,K.') 
disp('Input 0,0 for termination.') 
continue;
end

if (i==0)
stop_flag=1;
else   
% First we need to find out the user index in data
user_index = find(users==i);

% second we need to find out the movie index(s) that the user haven't seen
% (haven't rated.)

not_seen_in = find(rated_mat(user_index,:)==0);

p = zeros(length(not_seen_in),2);

for j=1:length(not_seen_in)
    movie_index = not_seen_in(j); % movie index
    p(j,1) = movies(movie_index); % store movie id for final report
    % Generate equation 1 of Breese et al (1998)
    Part1 = uv_bar(user_index);
    Part2 = (user_weight_mat(:,user_index)'*A1(:,movie_index))/...
    sum(abs(user_weight_mat(:,user_index))>0);
    
    % if the user voted for only one movie the denominator will be zero 
    % and the part 2 will be nan. In that case we use the same vote (part1) for all
    % movies because we don't have enough information about the user.
    
    if (isnan(Part2))
    Part2=0;
    end

    p(j,2)=Part1+Part2;
 
end

% return min of number of requested or number of available movie
report_num = min(length(not_seen_in),k);

% sort movies based on ranking
p = sortrows(p,-2);

% Print out the results
% Header
F1 = sprintf('%s%i','Recommended movies for user:',users(user_index));
disp(F1);
disp('----------------------------------------------');
for i=1:report_num

F2=sprintf('%s%s%i%s%2.2f',num2str(i),' - Movie ID : ', p(i,1), ' - Predicted Rate : ', p(i,2));
disp(F2);

end

disp('----------------------------------------------');
end

end












    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    