function [] = HMM_orientedGraph(type,gamma,prob,datasetName,pathFigures, HideZero, printOut, orientedOnly)
%FONCTION GRAPHIQUE Cette fonction permet d'afficher la structure HHMM
%   La fonction accepte les probabilités p1-p2-q1-q2
%   permettant de construire la matrice de transition du modèle, le nom du 
%   jeux de données pour donner le nom approprié au dossier enregistrer, 
%   ainsi que le chemin du dossier où enregistrer l'image
% ———————————————————————————————————————————
% Creation              : 21 janvier 2019
% ———————————————————————————————————————————
% Dernière modification : 25 février 2019
% ———————————————————————————————————————————
% Creation du nom du fichier a créer pour recevoir les résultats
newSubFolder = char(strcat(pathFigures,datasetName));

% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end


colors=[[1    0    0]; ...
        [0.9804,    0.5020,    0.4471]; ...
        [0.2000    0.6000    1.0000]; ...
        [0,         0,    0.6000]];

if strcmp(type,'HHMM') && (orientedOnly)
    
    gamma = round(gamma_Build_HHMM(prob),3);
    prob_round = round(prob,3);
    m = mean(mean(prob_round));

    links =[5 5 5 6 6 6 7 7 7 8 8 8 1 1 2 2 1 1 2 2 3 4 ; ...
            5 6 3 6 5 3 7 8 4 8 7 4 1 2 2 1 5 6 7 8 1 2];
    before_weights = [prob_round(1,:) prob_round(2,:) prob_round(3,:) prob_round(4,:) prob_round(5,1:2) prob_round(6,1:2) prob_round(7,1:2) prob_round(8,1:2) 1 1];
    
    %    1  2  3  4  5  6  7  8
    x = [1  6  1  6  0  2  5  7];
    y = [3  3  1  1  0  0  0  0];
    
    switch HideZero
        case 1
            s = [];
            t = [];
            weights = [];
            for i=1:length(before_weights)
                if before_weights(i)>0
                    s = [s links(1,i)];
                    t = [t links(2,i)];
                    weights = [weights before_weights(i)];
                end
            end

            G = digraph(s, t,weights,8);  

            n = length(G.Edges.Weight);
            edgeWeights = zeros(1,n);

            edgeWeights = strcat(string(100*G.Edges.Weight),'%');
            idxOut = findedge(G,[3 4],[1 2]);
            
            edgeWeights(idxOut) = ' ';
            
        case 0
            s = links(1,:);
            t = links(2,:);
            weights = before_weights

            G = digraph(s, t,weights,8);  

            n = length(G.Edges.Weight);
            edgeWeights = zeros(1,n);
            
            edgeWeights = strcat(string(100*G.Edges.Weight),'%');
            edgeWeights(9) = ' ';
            edgeWeights(10) = ' ';
    end

    % h=plot(G,'XData',x,'YData',y,'EdgeLabel', ...
    %     cellstr(strcat({'a1_11 = ', ...
    %             'a1_12 = ', ...
    %             'pi_1 = ', ...
    %             'pi_2 = ', ...
    %             'a1_21 = ', ...
    %             'a1_22 = ', ...
    %             'pi_3 = ', ...
    %             'pi_4 = ', ...
    %             '', ...
    %             '', ...
    %             'e1 = ', ...
    %             'a2_11 = ', ...
    %             'a2_12 = ', ...
    %             'e2 = ', ...
    %             'a2_21 = ', ...
    %             'a2_22 = ', ...
    %             'e3 = ', 'a2_33 = ', 'a2_34 = ', 'e4 = ', 'a2_43 = ', 'a2_44 = '},string(G.Edges.Weight).')),'LineWidth',2,'ArrowSize',15);

    h=plot(G,'XData',x,'YData',y,'EdgeLabel',cellstr(edgeWeights),'LineWidth',6*G.Edges.Weight/max(G.Edges.Weight)+0.0001,'ArrowSize',10);

    nRED = successors(G,5);
    nORANGE = successors(G,6);
    nGREEN = successors(G,7);
    nBLUE = successors(G,8);

    labelnode(h,[1 2 3 4 5 6 7 8],{'i_{1}' 'i_{2}' 'e_{1}' 'e_{2}' 's_{1}' 's_{2}' 's_{3}' 's_{4}'});
    nl = h.NodeLabel;
    h.NodeLabel = '';
    xd = get(h, 'XData');
    yd = get(h, 'YData');
    h.MarkerSize = 24;
    t = text(xd, yd, nl, 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Color',[0.99, 0.99 ,0.99]);
    t(6).Color = 'black';
    t(7).Color = 'black';

    highlight(h,5,'NodeColor',colors(1,:))
    highlight(h,6,'NodeColor',colors(2,:))
    highlight(h,7,'NodeColor',colors(3,:))
    highlight(h,8,'NodeColor',colors(4,:))
    highlight(h,5,nRED,'EdgeColor',colors(1,:))
    highlight(h,6,nORANGE,'EdgeColor',colors(2,:))
    highlight(h,7,nGREEN,'EdgeColor',colors(3,:))
    highlight(h,8,nBLUE,'EdgeColor',colors(4,:))
    
    % Retrait de axes
    set(gca,'XTick',[],'YTick',[])
%     set(gca,'visible','off')

    if printOut
        title(strcat('Graphique Hiérarchie : Modèle ',{' '},type,{' '},'pour',{' '},datasetName))
        set(gcf,'PaperPositionMode','auto',...
            'Units', 'Normalized',...
            'OuterPosition', [0, 0.04, 1, 0.96])

        print(strcat(newSubFolder,'\',datasetName,'_Oriented_Graph_HHMM.png'),'-dpng','-r500')

        close(gcf)
    end
end

%  ———————————————————————————————————————————
%  ———————————————————————————————————————————
if ~orientedOnly
    links = [];
    prob = [];
    nbRegime = length(gamma);
    gamma =round(gamma,3);

    for i=1:nbRegime
        for j=1:nbRegime
            if gamma(i,j)>0
                links = [links [i ; j]];
                prob = [prob gamma(i,j)];
            end
        end
    end

    %    1  2  3  4  5  6  7  8
    x = [0 0 3 3];
    y = [3 0 0 3];

    s = links(1,:);
    t = links(2,:);
    weights = prob;

    G = digraph(s, t,weights,4);  
    % 
    % % h=plot(G,'XData',x,'YData',y,'EdgeLabel', ...
    % %     cellstr(strcat({'a1_11 = ', ...
    % %             'a1_12 = ', ...
    % %             'pi_1 = ', ...
    % %             'pi_2 = ', ...
    % %             'a1_21 = ', ...
    % %             'a1_22 = ', ...
    % %             'pi_3 = ', ...
    % %             'pi_4 = ', ...
    % %             '', ...
    % %             '', ...
    % %             'e1 = ', ...
    % %             'a2_11 = ', ...
    % %             'a2_12 = ', ...
    % %             'e2 = ', ...
    % %             'a2_21 = ', ...
    % %             'a2_22 = ', ...
    % %             'e3 = ', 'a2_33 = ', 'a2_34 = ', 'e4 = ', 'a2_43 = ', 'a2_44 = '},string(G.Edges.Weight).')),'LineWidth',2,'ArrowSize',15);
    % n = length(G.Edges.Weight);
    % edgeWeights = zeros(1,n);
    % 
    % edgeWeights = string(G.Edges.Weight);
    % edgeWeights(9) = '';
    % edgeWeights(10) = '';
    % 
    h=plot(G,'XData',x,'YData',y,'EdgeLabel',cellstr(strcat(string(100*G.Edges.Weight),'%')),'LineWidth',7*G.Edges.Weight/max(G.Edges.Weight)+0.0001,'ArrowSize',10);


    nRED = successors(G,1);
    nORANGE = successors(G,2);
    nGREEN = successors(G,3);
    nBLUE = successors(G,4);


    labelnode(h,[1 2 3 4],{'s_{1}' 's_{2}' 's_{3}' 's_{4}'});
    nl = h.NodeLabel;
    h.NodeLabel = '';
    xd = get(h, 'XData');
    yd = get(h, 'YData');
    h.MarkerSize = 24;
    t = text(xd, yd, nl, 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Color',[0.99, 0.99 ,0.99]);
    t(2).Color = 'black';
    t(3).Color = 'black';

    highlight(h,1,'NodeColor',[1    0    0])
    highlight(h,2,'NodeColor',[0.9804,    0.5020,    0.4471])
    highlight(h,3,'NodeColor',[0,    1,    1])
    highlight(h,4,'NodeColor',[0, 0, 1])
    highlight(h,1,nRED,'EdgeColor',[1    0    0])
    highlight(h,2,nORANGE,'EdgeColor',[0.9804,    0.5020,    0.4471])
    highlight(h,3,nGREEN,'EdgeColor',[0,    1,    1])
    highlight(h,4,nBLUE,'EdgeColor',[0, 0, 1])

    % Retrait de axes
    set(gca,'XTick',[],'YTick',[])
    % set(gca,'visible','off')

    if printOut
        title(strcat('Processus Changement de Régime : Modèle ',{' '},type,{' '},'pour',{' '},datasetName))
        set(gcf,'PaperPositionMode','auto',...
            'Units', 'Normalized',...
            'OuterPosition', [0, 0.04, 1, 0.96])

        print(strcat(newSubFolder,'\',datasetName,'_Regime_Switching_Graph_',type,'.png'),'-dpng','-r500')

        close(gcf)
    end
end

end

