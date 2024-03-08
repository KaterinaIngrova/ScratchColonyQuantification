%% Studium vlivu ropivakain hydrochloridu na migrační a proliferační schopnost nádorových buněk lidského osteosarkomu SaOS-2
% Kateřina Ingrová - BTBIO, 219244
% MATLAB R2020a

%HODNOCENÍ VLIVU LA ROPIVAKAINU NA MIGRACI BUNĚK SaOS-2

results = Scratch_function('C:\Users\user\Desktop\diplomka\matlab\nove bunky snimky Inna\Scratch', 0.6369);

function [results] = Scratch_function(cesta_k_obrazkum, velikost_pixelu_mikrom);
% 1. Získání seznamu obrázků
    cd(cesta_k_obrazkum);
    images = dir(fullfile(cesta_k_obrazkum, '*.tif'));

    % 2. Inicializace výsledků
    number_of_images = length(images);
    results = cell(number_of_images, 3);
    results(1,1) = {'Název obrázku'};
    results(1,2) = {'Průměrná šířka rýhy [μm]'};
    results(1,3) = {'Plocha rýhy [mm2]'};

    % 3. Zpracování jednotlivých obrázků
    %For cyklus k procházení jednotlivých snímků.
    for i=9:9 %number_of_images 
        name = images(i).name; %název snímku
        results(i+1,1)={name}; %uložení názvu do pole
        concrete_image = imread(name); 
        disp('------------------------');
        disp(['Název snímku ', name]);

        %Zobrazení snímku konkretního snímku ze složky.
        figure(1); 
        subplot(2,3,1)
        imshow(concrete_image);
        title('1');

        %Zjištění rozměrů snímku.
        size_of_image = size(concrete_image); 
        number_of_rows = size_of_image(1);
        number_of_columns = size_of_image(2);
        number_of_pixels = number_of_rows*number_of_columns;

        %Převednení šedotónového snímku do double formátu.
        double_gray = double(concrete_image);

        %Zlepšení kontrastu snímku.
        contrast_enhancement = imadjust(concrete_image);

        %Nalezení hran pomocí Cannyho detektoru.
        edge_image = edge(contrast_enhancement, 'Canny', [0.08 0.1], 5); %0,03

        %Zobrazení hranové reprezentace.
        subplot(2,3,2)
        imshow(edge_image);
        title('2');

        %Nalezení kruhů, které odpovídají mrtvým buňkám.
        [centers, radii, metric] = imfindcircles(edge_image,[10 25]);

        %Zobrazení mrvých buněk.
%         figure(4);
%         imshow(edge_image); 
%         viscircles(centers, radii, 'EdgeColor', 'b'); %zobrazení kruhů modře
        outputImage = edge_image; 

        %Vytvoření jasového profilu snímku.
        imageProfile = mean(double_gray, 1);  % Průměr hodnot pro každý sloupec. contrast_enhanement

        %Zobrazení jasového profilu.
%         figure(5);
%         plot(imageProfile, 'LineWidth', 2);
%         title('Jasový profil pro každý sloupec');
%         xlabel('Pozice');
%         ylabel('Průměrný jas');

        %Nalezení pozice minima v jasovém profilu.
        [minValue, minIndex] = min(imageProfile);

        %Zvýraznění minima na jasovém profilu.
        %hold on;
        %plot(minIndex, minValue, 'ro', 'MarkerSize', 10);  % Červený kruh označuje minimum.


        %fprintf('Minimum jasového profilu je na pozici: %d\n', minIndex);

        %Filtrace maxim s průměrným jasem.
        [pks, locs] = findpeaks(imageProfile, 'MinPeakHeight', mean(imageProfile)); %-3
        maxValues = pks;
        maxIndices = locs;

        %Nalezení nejbližšího maxima vlevo od minima
        maxima_vlevo = locs(locs < minIndex);
        nearest_maximum_left = max(maxima_vlevo);

        %Nalezení nejbližšího maxima vpravo od minima
        maxima_vpravo = locs(locs > minIndex);
        nearest_maximum_right = min(maxima_vpravo);

        % Výsledky
        %disp(['Nejbližší maximum vlevo: ', num2str(nearest_maximum_left)]);
        %disp(['Nejbližší maximum vpravo: ', num2str(nearest_maximum_right)]);

        %plot(nearest_maximum_left, imageProfile(nearest_maximum_left), 'go', 'MarkerSize', 8);  % Zelený kruh označuje maximum.
        %plot(nearest_maximum_right, imageProfile(nearest_maximum_right), 'go', 'MarkerSize', 8);  % Zelený kruh označuje maximum.

        %hold off


        %Hranice průměrné hodnoty jasu pro odstranění buňek (upravte podle potřeby).
        brightnessThreshold = 70; %70

        %Odstranění mrtvých buněk z oblasti hranice rýhy.
        for h = 1:size(centers,1) 
            center = round(centers(h, :));
            radius = round(radii(h)); %4

            %Vytvoření masky pro daný kruh.
            [xx, yy] = meshgrid(1:size(outputImage, 2), 1:size(outputImage, 1));
            mask = (xx - center(1)).^2 + (yy - center(2)).^2 <= radius^2;

            %Omezení masky na sloupce od 640 do 1280 - na rýhu.
            mask(:, 1:nearest_maximum_left) = false; %768
            mask(:, nearest_maximum_right:end) = false; %1152

            % Průměrná hodnota jasu pro buňky v masce.
            averageBrightness = mean(contrast_enhancement(mask));

            %Odstranění mrtvých buněk z binární hranové reprezentace.
            % Odstranění buňek s průměrnou hodnotou jasu nad hranicí.
            if averageBrightness > brightnessThreshold %%|| any(surroundingMask(:))
                outputImage(mask) = 0; % Nastavení hodnoty na 0 (černá) pro mrtvé buňky.
            end
        end

        %Zobrazení snímků po odstranění mrtvých buněk z oblasti rýhy.
        subplot(2,3,3)
        imshow(outputImage);
        title('3');

        %Použití morfologických operátorů. 
        se=strel('disk',18); %20 puvodne, pozdeji 18
        after=imclose(outputImage,se); %aplikace closingu
        se1 = strel('disk',8);%5 puvodne
        after2=imopen(after,se1); %aplikace openingu

        %Zobrazení snímku po použítí morfologických operací.
        subplot(2,3,4)
        imshow(after2);
        title('4');
        
        copy_2 = after2;

        %Nastavení pro prvních 5 a posledních 5 řádků první a poslední třetiny
        %sloupců na 1. 
    %     for row = 1:number_of_rows
    %         for column = 1:number_of_columns
    %             if (column <= number_of_columns / 3 || column > (2 * number_of_columns / 3)) && ((row <= 5 || row > (number_of_rows - 5)) || (column <= 5 || column > (number_of_columns - 5)))
    %                 copy_2(row, column) = 1;
    %             end
    %         end
    %     end

        for row = 1:number_of_rows
            for column = 1:number_of_columns
                if (column <=  nearest_maximum_left || column > nearest_maximum_right) && ((row <= 5 || row > (number_of_rows - 5)) || (column <= 5 || column > (number_of_columns - 5)))
                    copy_2(row, column) = 1;
                end
            end
        end

        number_white = sum(copy_2(:) == 1);
        number_black = sum(copy_2(:) == 0);

        %Vyplnění děr na vysegmentovaném obrázku rýhy.
        segmented_image = imfill(copy_2,'holes');

        %Zabránění vyplňování děr v oblasti rýhy, pokud poměr černých pixelů
        %bude menší než 10 %. 
        if number_black/number_of_pixels < 0.1 %20
            %start_column = nearest_maximum_left; %number_of_columns / 3;
            %end_column = nearest_maximum_right; %2 *(number_of_columns / 3);
            %segmented_image(:, 670:1300) = copy_2(:, 670:1300);

            % lokalni mapa 
            lokalni_std = stdfilt(contrast_enhancement);

            %Zobrazení lokální mapy směrodatných odchylek.
%             figure(8);
%             imshow(lokalni_std, []);

            local_entropy = entropyfilt(lokalni_std);

            %Zobrazení lokální mapy entropie.
%             figure(9);
%             imshow(local_entropy, []);

            image_2 = zeros(size(local_entropy));

            % Nastavení hodnot pixelů na 1, pokud lokální entropie je menší než 1
            image_2(local_entropy > 0) = 1;
            image_2 = imfill(image_2,'holes');
            image_2 = ~image_2;

            %Zobrazení binárního obrazu
%             figure(10);
%             imshow(image_2);


            se=strel('disk',12); %20 puvodne, pozdeji 18
            o=imclose(image_2,se); %aplikace closingu
            se1 = strel('disk',3);%5 puvodne
            after_o=imopen(o,se1); %aplikace openingu

            %Zobrazení obrazu po morfologických operacích.
            subplot(2,3,4)
            imshow(after_o);
            title('4');

            for row = 1:number_of_rows
                for column = 1:number_of_columns
                    if (column <=  nearest_maximum_left || column > nearest_maximum_right) && ((row <= 5 || row > (number_of_rows - 5)) || (column <= 5 || column > (number_of_columns - 5)))
                        after_o(row, column) = 1;
                    end
                end
            end

            segmented_image = imfill(after_o, 'holes');

            segmented_image(:, 670:1300) = after_o(:, 670:1300);
            
            % Najdi oblasti s hodnotou nula v segmented_image
            zero_regions = bwconncomp(segmented_image == 0);

            % Zobraz jednotlivé oblasti s nulovými hodnotami
%             figure;
%             imshow(segmented_image);
%             hold on;

            for i = 1:zero_regions.NumObjects
                % Získání pixelů pro každou oblast
                pixels = zero_regions.PixelIdxList{i};

                % Převedení lineárního indexu na souřadnice
                [rows, cols] = ind2sub(size(segmented_image), pixels);

                % Vykresli hranice oblasti
%                 plot(cols, rows, 'r.', 'MarkerSize', 5);
            end
% 
%             hold off;
%             title('Oblasti s nulovými hodnotami');
            
            % Vytvoř pole pro ukládání průměrných hodnot kontrastu v jednotlivých oblastech
            average_contrast_values = zeros(1, zero_regions.NumObjects);

            for i = 1:zero_regions.NumObjects
                % Získání pixelů pro každou oblast
                pixels = zero_regions.PixelIdxList{i};

                % Extrahuj hodnoty kontrastu v dané oblasti
                contrast_values = contrast_enhancement(pixels);

                % Vypočti průměrnou hodnotu kontrastu v dané oblasti
                average_contrast = mean(contrast_values);

                 % Pokud je průměrná hodnota kontrastu větší než threshold
                if average_contrast > 70
                    % Nastav bílou barvu v segmented_image pro danou oblast
                    segmented_image(pixels) = 1;  % 1 odpovídá bílé barvě v binárním obrázku
                end
            end
        end



        %Zobrazení vysegmentovaného obrázku.
        subplot(2,3,5)
        imshow(segmented_image);
        title('5');

        %Výpočet šířky rýhy kromě prvních a posledních 5 řádků.
        width_scratch_each_row = sum(segmented_image(6:(end-5),:) == 0, 2);
        mean_of_scratch = mean(width_scratch_each_row);
        mean_of_scratch_mikrometr = mean_of_scratch*velikost_pixelu_mikrom; % v mm (1447 1918 - odpovídá 300 mikrometrum což je 471 => 1 pixel = 0,6369)
        %sirka_ryhy = zeros(number_of_rows, number_of_columns);
        %sirka = number_of_columns/3;


    %     for radky = 1:number_of_rows;
    %         for column = sirka:(sirka*2);
    %             if segmented_image(radky,column) ==0;
    %                sirka_ryhy(radky,column)= 1;
    %             end
    %         end
    %     end
    %     
    %     
        %suma_ryhy = sum(sirka_ryhy,2);
        %rows_to_remove = any(suma_ryhy == 0,2); %%odstranit prvnich a poslednich 5 najit indexy
        %suma_ryhy(rows_to_remove, :) = []; %odstranit
        %prumer_ryhy = mean(suma_ryhy);  

        %Uložení průměrné rýhy v mikrometech do pole.
        results(i+1,2)={mean_of_scratch_mikrometr};

        edgeImage = bwperim(segmented_image);

        % Vyloučit první a poslední řádek a sloupec
        edgeImage([1:5, number_of_rows-4:number_of_rows], :) = 0;
        edgeImage(:, [1, end]) = 0;
        % Získat pozice bílých pixelů na hranovém obrazu
        [row, col] = find(edgeImage == 1);

        % Přidat hranu rýhy na původní obrazek
        subplot(2,3,6)
        imshow(concrete_image);
        title('6');
        hold on;
        plot(col, row, 'g.', 'Marker', '.', 'MarkerSize', 5); % Jasné modré body
        hold off;
        %title('Hranice vysegmentované rýhy na původním obrazku');

        disp(['Průměrná šířka rýhy: ', num2str(mean_of_scratch_mikrometr), ' μm']);


        % Výpočet plochy rýhy
        plocha_ryhy = (sum(width_scratch_each_row)*(velikost_pixelu_mikrom^2)/1000000); % Předpokládáme, že šířka rýhy je součtem šířek v každém řádku.

        % Případně, pokud máte šířku rýhy v každém řádku ve vektoru sirka_rhy_v_kazdem_radku, můžete použít:
        % plocha_rhy = sum(sirka_rhy_v_kazdem_radku) * pocet_sloupcu;

        disp(['Plocha rýhy: ' num2str(plocha_ryhy) ' mm^2']);
        results(i+1,3)={plocha_ryhy};

    end
end

