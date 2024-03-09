%The effect of gentamicin on the migration of SaOS-2 cell line

%Example of calling function with the obtained results
results = Scratch_function('C:\Users\user\Desktop\diplomka\matlab\nove bunky snimky Inna\Scratch', 0.6369);

%The function for the evaluation - the inputs are a path to the microscopy images and size of
%pixel in micrometers, the obtained results are saved into cell (scratch
%width in micrometers and the scratch area in mm2)

function [results] = Scratch_function(path_to_image, size_of_pixel_micrometer);
    
    %List of images 
    cd(path_to_image);
    files = dir(fullfile(path_to_image, '*.*'));
    images = files(~[files.isdir]);

    %Inicialization of cell for storing results
    number_of_images = length(images);
    results = cell(number_of_images, 3);
    results(1,1) = {'Name of the image'};
    results(1,2) = {'Average scratch width [μm]'};
    results(1,3) = {'Scratch area [mm2]'};

    %For loop for loading current image
    %For cyklus k procházení jednotlivých snímků.
    for i=1:number_of_images
        name = images(i).name; %name of image
        results(i+1,1)={name}; %saving the name of image to created cell
        concrete_image = imread(name); %load concrete image
        disp('------------------------');
        disp(['Název snímku ', name]);

        %Suplot of important steps of algorithm 
        figure(1); 
        subplot(2,3,1)
        imshow(concrete_image); 
        title('1');

        %Determination of image size
        size_of_image = size(concrete_image); 
        number_of_rows = size_of_image(1);
        number_of_columns = size_of_image(2);
        number_of_pixels = number_of_rows*number_of_columns;
        
        %Finding the dimension of image, in the case of RGB image is image
        %converted to grayscale, in both cases is applying contrast
        %enhancement on images
        if ndims(concrete_image) == 3 
        sedoton = rgb2gray(concrete_image);
        contrast_enhancement = imadjust(sedoton);
    
        elseif ndims(concrete_image) == 2
        contrast_enhancement = imadjust(concrete_image);
        
        end

        
        %Finding the edges of image by Canny detector
        edge_image = edge(contrast_enhancement, 'Canny', [0.08 0.1], 5); 

        %Suplot of important steps of algorithm 
        subplot(2,3,2)
        imshow(edge_image);
        title('2');

        %Find circles with coresponding to dead cells
        [centers, radii, metric] = imfindcircles(edge_image,[10 25]);
        % figure(4);
        % imshow(edge_image); 
        % viscircles(centers, radii, 'EdgeColor', 'b'); 
        
        %Make copy of edge image for segmentation
        outputImage = edge_image; 

        %Make brightness profile of image
        imageProfile = mean(contrast_enhancement, 1);
        % figure(5);
        % plot(imageProfile, 'LineWidth', 2);
        % title('Brightness profile');
        % xlabel('Pozice');
        % ylabel('Průměrný jas');

        %Finding minimum in brightness profile
        [minValue, minIndex] = min(imageProfile);
        %hold on;
        %plot(minIndex, minValue, 'ro', 'MarkerSize', 10);

        %Filtration of maximas with average brightness
        [pks, locs] = findpeaks(imageProfile, 'MinPeakHeight', mean(imageProfile));
        maxValues = pks;
        maxIndices = locs;

        %Finding the nearest maxima on the left from minima, witch represent
        %a border of scratch
        left_maxima = locs(locs < minIndex);
        nearest_maximum_left = max(left_maxima);

        %Finding the nearest maxima on the right from minima, witch
        %represent a border of scratch
        right_maxima = locs(locs > minIndex);
        nearest_maximum_right = min(right_maxima);
        %plot(nearest_maximum_left, imageProfile(nearest_maximum_left), 'go', 'MarkerSize', 8);  
        %plot(nearest_maximum_right, imageProfile(nearest_maximum_right), 'go', 'MarkerSize', 8); 
        %hold off

        %Threshold of average brightness of dead cells to removing them
        brightnessThreshold = 70; 

        %Removing dead cells from the scratch region
        for h = 1:size(centers,1) 
            center = round(centers(h, :));
            radius = round(radii(h)); 

            %Create mask for each dead cell
            [xx, yy] = meshgrid(1:size(outputImage, 2), 1:size(outputImage, 1));
            mask = (xx - center(1)).^2 + (yy - center(2)).^2 <= radius^2;

            %Remove them only from the scratch region
            mask(:, 1:nearest_maximum_left) = false; 
            mask(:, nearest_maximum_right:end) = false; 

            %Average brightness of dead cells
            averageBrightness = mean(contrast_enhancement(mask));

            %Removing dead cells from edge representation, the dead cells
            %must have average brightness higher then threshold brightness
            if averageBrightness > brightnessThreshold 
                outputImage(mask) = 0; %set death cells to zero
            end
        end

        %Suplot of important steps of algorithm 
        subplot(2,3,3)
        imshow(outputImage);
        title('3');

        %Using morphological operations for segmentation 
        se=strel('disk',18); 
        after=imclose(outputImage,se);
        se1 = strel('disk',8);
        after2=imopen(after,se1);
        
        %Suplot of important steps of algorithm 
        subplot(2,3,4)
        imshow(after2);
        title('4');
        
        
        morphological_operations = after2;

        %Setting first and last 5 rows to 1 for better hole filling 
        for row = 1:number_of_rows
            for column = 1:number_of_columns
                if (column <=  nearest_maximum_left || column > nearest_maximum_right) && ((row <= 5 || row > (number_of_rows - 5)) || (column <= 5 || column > (number_of_columns - 5)))
                    morphological_operations(row, column) = 1;
                end
            end
        end
        
        %Calculate numbef of pixels witch represent cells and pixels witch
        %represent scratch
        number_white = sum(morphological_operations(:) == 1);
        number_black = sum(morphological_operations(:) == 0);

        %Hole filling
        segmented_image = imfill(morphological_operations,'holes');

        %Diffent way (hole filling was problem) for images wich have got
        %less then 10% of sctach pixels
        if number_black/number_of_pixels < 0.1 

            %Create local map of standard deviation
            local_std = stdfilt(contrast_enhancement);
            %figure(8);
            %imshow(local_std, []);

            %Create local map of entropy
            local_entropy = entropyfilt(local_std);
            %figure(9);
            %imshow(local_entropy, []);

            %Create zeros mask for segmentation
            segmentation_second_case = zeros(size(local_entropy));

            %First part of segmenation - put to zero mask 1 if
            %local_entropy is higher than 0
            segmentation_second_case(local_entropy > 0) = 1;
            segmentation_second_case = imfill(segmentation_second_case,'holes');
            segmentation_second_case = ~segmentation_second_case;
            %figure(10);
            %imshow(segmentation_second_case);

            %Another part - use morphological operations
            se=strel('disk',12); 
            o=imclose(segmentation_second_case,se); 
            se1 = strel('disk',3);
            after_o=imopen(o,se1); 

            %Suplot of important steps of algorithm 
            subplot(2,3,4)
            imshow(after_o);
            title('4');
            
            %Put first 5 and last 5 rows to 1 for better segmentation
            for row = 1:number_of_rows
                for column = 1:number_of_columns
                    if (column <=  nearest_maximum_left || column > nearest_maximum_right) && ((row <= 5 || row > (number_of_rows - 5)) || (column <= 5 || column > (number_of_columns - 5)))
                        after_o(row, column) = 1;
                    end
                end
            end
            
            %Fill holes 
            segmented_image = imfill(after_o, 'holes');
            
            %Scratch area without filling holes
            segmented_image(:, nearest_maximum_left:nearest_maximum_right) = after_o(:, nearest_maximum_left:nearest_maximum_right);
        
            
            %Find in segmented image scratch part
            zero_regions = bwconncomp(segmented_image == 0);
            %figure;
            %imshow(segmented_image);
            %hold on;
            
            %For loop for all sratch part
            for i = 1:zero_regions.NumObjects
                %Get pixels of each part
                pixels = zero_regions.PixelIdxList{i};
                [rows, cols] = ind2sub(size(segmented_image), pixels);
                %plot(cols, rows, 'r.', 'MarkerSize', 5);
            end
            %hold off;
            %title('Oblasti s nulovými hodnotami');
            
            %Find average brightness for them
            average_brightness_value = zeros(1, zero_regions.NumObjects);

            for i = 1:zero_regions.NumObjects
                %Find pixels for each region
                pixels = zero_regions.PixelIdxList{i};

                %Find brightness in region
                brightness_values = contrast_enhancement(pixels);

                %Calcuated the mean value in region
                average_brightness = mean(brightness_values);

                %If the average brightness is higher then threshold, it is
                %dead cell and fill it 
                if average_brightness > 70
                    segmented_image(pixels) = 1; 
                end
            end
        end



        %Suplot of important steps of algorithm 
        subplot(2,3,5)
        imshow(segmented_image);
        title('5');

        %Calculating the parameters of scratch
        width_scratch_each_row = sum(segmented_image(6:(end-5),:) == 0, 2);
        mean_of_scratch = mean(width_scratch_each_row);
        mean_of_scratch_micrometer = mean_of_scratch*size_of_pixel_micrometer; 
        area_scratch = (sum(width_scratch_each_row)*(size_of_pixel_micrometer^2)/1000000); % Předpokládáme, že šířka rýhy je součtem šířek v každém řádku.
        disp(['Average scratch width: ', num2str(mean_of_scratch_micrometer), ' μm']);
        disp(['Scratch area: ' num2str(area_scratch) ' mm^2']);

        %Saving the results 
        results(i+1,2)={mean_of_scratch_micrometer};
        results(i+1,3)={area_scratch};

        %Find the boundaries of scratch on original image
        edgeImage = bwperim(segmented_image);
        edgeImage([1:5, number_of_rows-4:number_of_rows], :) = 0;
        edgeImage(:, [1, end]) = 0;
        [row, col] = find(edgeImage == 1);

        %Suplot of important steps of algorithm 
        subplot(2,3,6)
        imshow(concrete_image);
        title('6');
        hold on;
        plot(col, row, 'g.', 'Marker', '.', 'MarkerSize', 5);
        hold off;

    end
end


