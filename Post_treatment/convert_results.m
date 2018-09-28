filename='FINAL_error_IP__p3';
output=strcat(filename,'_output.txt');
file=fopen(output,'w');
result=struct2array(load(strcat(filename,'.mat')));
[N1,N2]=size(result);
fprintf(file,'\\begin{table}[H]\\centering\\makebox[\\textwidth][c]{\\begin{tabular}{');
for j=1:1:N2
    if j==1
        fprintf(file,'|c');
    elseif j==N2
        fprintf(file,'|c|}');
    else
        fprintf(file,'|c');
    end
end
N3=4;
for i=1:1:N1
    for j=1:1:N2
        switch j
            case 1
                switch i
                    case 1
                        fprintf(file,'%d &',result(i,j));
                    case 2
                        fprintf(file,'%4.3f &',result(i,j));
                    otherwise
                        quotient=floor((i-3)/N3);
                        res=(i-3)-N3*quotient;
                        switch quotient
                            case 0
                            switch res
                                case 0
                                    fprintf(file,'\\multirow{4}{*}{$c_{DG}$} &');
                                case 1
                                    fprintf(file,'&');
                                case 2
                                    fprintf(file,'&');
                                case 3
                                    fprintf(file,'&');                                               
                            end
                            case 1
                            switch res
                                case 0
                                    fprintf(file,'\\multirow{4}{*}{$c_{SD}$} &');
                                case 1
                                    fprintf(file,'&');
                                case 2
                                    fprintf(file,'&');
                                case 3
                                    fprintf(file,'&');                                               
                            end
                            case 2
                            switch res
                                case 0
                                    fprintf(file,'\\multirow{4}{*}{$c_{HU}$} &');
                                case 1
                                    fprintf(file,'&');
                                case 2
                                    fprintf(file,'&');
                                case 3
                                    fprintf(file,'&');                                               
                            end
                            case 3
                            switch res
                                case 0
                                    fprintf(file,'\\multirow{4}{*}{$c_{+}$} &');
                                case 1
                                    fprintf(file,'&');
                                case 2
                                    fprintf(file,'&');
                                case 3
                                    fprintf(file,'&');                                               
                            end
                        end
                end
            case 2
                switch i
                    case 1
                        
                    case 2
                        fprintf(file,'%4.3f &',result(i,j));
                    otherwise
                        quotient=floor((i-3)/N3);
                        res=i-3-N3*quotient
                        switch res
                            case 0
                                fprintf(file,'$\\kappa_{DG}$&');
                            case 1
                                fprintf(file,'$\\kappa_{SD}$&');
                            case 2
                                fprintf(file,'$\\kappa_{HU}$&');
                            case 3
                                fprintf(file,'$\\kappa_{+}$&');                                               
                        end
                end
            case N2
                switch i
                    case 1
                        fprintf(file,'%d \\\\  \\hline ',result(i,j));
                    case 2
                        fprintf(file,'%4.3f \\\\  \\hline ',result(i,j));
                    case N1
                        fprintf(file,'%8.2e \\\\ \\hline ',result(i,j));
                    otherwise
                        fprintf(file,'%8.2e \\\\  \\hline ',result(i,j));
                end                
                
            otherwise
                switch i
                    case 1
                        fprintf(file,'%d &',result(i,j));
                    case 2
                        fprintf(file,'%4.3f &',result(i,j));
                    otherwise
                        fprintf(file,regexprep(sprintf('%8.2e &',result(i,j)),'e-','e-'));
                end                
                
        end
    end
end
fprintf(file,'\\end{tabular}\n}\n\\caption{$s_{numerical}^{*}$ (BR2) for $N_{p}=3$ for 32 elements}\n\\label{result BR2 Np=3}\n\\end{table}');


