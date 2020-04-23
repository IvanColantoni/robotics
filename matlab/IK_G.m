%metodo del Gradiente a passo costante per soluzione di problemi di
%cinematica inversa. E' necessario settare il passo dell'iterata per il
%metodo del gradiente "alpha", fornire in input il vettore symbolico
%P(q1,...,qn), il vettore f_c della configurazone finale dell'end effector
% e una variabile char se si vuole visualizzare l'andamento
%della convergenza dei due metodi 
function[]=IK_G(pq,q10,q20,alpha,f_c,varargin)
 syms q1 q2 ;
 Jx=(jacobian(pq)).';
 %Q=[q1,q2];
 q1=q10;
 q2=q20;
 e1=(f_c-eval(pq));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
 %disp(e1)
 fileID= fopen('IK_G.txt', 'w');
 fprintf(fileID,'%12s %12s %12s %12s\n','q1_G','q2_G');%, 'q1_G','q2_G');
 NORM=1:100;
 figure(1);
 for k=1:100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
   % disp(eval(Jx));
    G=(eval(Jx))*e1.';
    %disp(alpha*G);
    q1=q1+alpha*G(1);
    q2=q2+alpha*G(2);
    e1=f_c-eval(pq);
    NORM(k)=norm(e1,2);
    fprintf(fileID, '\n%12.5f %12.5f\n', q1, q2);%, x1, x2);
  end
plot(1:100,NORM,'b'); hold on
fclose(fileID);

