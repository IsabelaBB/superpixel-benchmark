/* -- Serge Bobbia : serge.bobbia@u-bourgogne.fr -- Le2i 2018
 * This work is distributed for non commercial use only,
 * it implements the IBIS method as described in the CVPM 2018 paper.
 *
 * Read the ibis.h file for options and benchmark instructions
 */

#include "ibis.h"

IBIS::IBIS(int _maxSPNum, int _compacity ) {
    labels = nullptr;
    maxSPNumber = _maxSPNum;
    compacity = _compacity;
    size = 0;

    // memory allocation
    Xseeds = new float[maxSPNumber];
    Yseeds = new float[maxSPNumber];
    lseeds = new float[maxSPNumber];
    aseeds = new float[maxSPNumber];
    bseeds = new float[maxSPNumber];

    Xseeds_init = new float[maxSPNumber];
    Yseeds_init = new float[maxSPNumber];

    Xseeds_Sum = new float[maxSPNumber];
    Yseeds_Sum = new float[maxSPNumber];
    lseeds_Sum = new float[maxSPNumber];
    aseeds_Sum = new float[maxSPNumber];
    bseeds_Sum = new float[maxSPNumber];

    countPx = new float[maxSPNumber];

    inv = new float[maxSPNumber];
    adjacent_sp = new int[size_roi*maxSPNumber];
    count_adjacent = new int[maxSPNumber];

}

IBIS::~IBIS() {
    delete[] labels;
    delete[] avec;
    delete[] lvec;
    delete[] bvec;
    delete[] inv;

    // allocated in constructor
    delete[] Xseeds;
    delete[] Yseeds;
    delete[] lseeds;
    delete[] aseeds;
    delete[] bseeds;

    delete[] Xseeds_init;
    delete[] Yseeds_init;

    delete[] Xseeds_Sum;
    delete[] Yseeds_Sum;
    delete[] lseeds_Sum;
    delete[] aseeds_Sum;
    delete[] bseeds_Sum;

    delete[] countPx;

    delete[] mask_buffer;
    delete[] mask_size;
    delete[] x_vec;
    delete[] y_vec;
    delete[] vertical_index;

    delete[] adjacent_sp;
    delete[] count_adjacent;
    delete[] initial_repartition;
    delete[] processed;

}

void IBIS::initSeeds() {
    int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip, yerrperstrip;
    int xoff, yoff;
    int x, y;
    int xe, ye, xe_1, ye_1;
    int start_y, final_y, start_x, final_x;
    int* tmp_adjacent_sp = new int[size_roi*maxSPNumber];
    int* tmp_count_adjacent = new int[maxSPNumber];

    xstrips = width / SPTypicalLength;
    ystrips = height / SPTypicalLength;

    xerr = width - SPTypicalLength * xstrips;
    yerr = height - SPTypicalLength * ystrips;

    xerrperstrip = (double)xerr / xstrips;
    yerrperstrip = (double)yerr / ystrips;

    xoff = SPTypicalLength / 2;
    yoff = SPTypicalLength / 2;

    n = 0;
    for (y = 0; y < ystrips; y++)
    {
        ye = (int)(y*yerrperstrip);
        ye_1 = (int)((y+1)*yerrperstrip);

        int seedy = y * SPTypicalLength + yoff + ye;

        if( y == 0 ) {
            start_y = 0;
            final_y = SPTypicalLength + ye_1;

        }
        else {
            start_y = y * SPTypicalLength + ye;
            final_y = ( (y + 1) * SPTypicalLength + ye_1 >= height ) ? height-1 : (y + 1) * SPTypicalLength + ye_1;

        }

        for (x = 0; x < xstrips; x++)
        {
            int seedx;
            xe = (int)(x*xerrperstrip);
            xe_1 = (int)((x+1)*xerrperstrip);
            seedx = x * SPTypicalLength + xoff + xe;

            if( x == 0 ) {
                start_x = 0;
                final_x = SPTypicalLength + xe_1;

            }
            else {
                start_x = x * SPTypicalLength + xe;

                final_x = ( (x + 1) * SPTypicalLength + xe_1 > width ) ? width : (x + 1) * SPTypicalLength + xe_1;

            }

            Xseeds_init[n] = (float) seedx;
            Yseeds_init[n] = (float) seedy;

            // fill line by line
            for( int index_y=start_y; index_y<=final_y; index_y++ ) {
                std::fill( initial_repartition + index_y*width + start_x, initial_repartition + index_y*width + final_x, n );

            }

            // list adjacents seeds
            tmp_count_adjacent[n] = 0;
            for( int roi_y=-(sqrt(size_roi) - 1)/2; roi_y <= (sqrt(size_roi) - 1)/2; roi_y++ ) {
                for( int roi_x=-(sqrt(size_roi) - 1)/2; roi_x <= (sqrt(size_roi) - 1)/2; roi_x++ ) {
                    if( !( y + roi_y < 0 || y + roi_y >= ystrips || x + roi_x < 0 || x + roi_x >= xstrips ) ) {
                        tmp_adjacent_sp[size_roi*n+tmp_count_adjacent[n]] = n + roi_y*xstrips + roi_x;
                        tmp_count_adjacent[n]++;

                    }

                }

            }

            n++;
        }
    }
    SPNumber = n;

    for(int i=0; i<SPNumber; i++) {
        count_adjacent[i] = 0;

        for(int j=0; j<tmp_count_adjacent[i]; j++) {
            if( tmp_adjacent_sp[size_roi*i+j] >= 0 && tmp_adjacent_sp[size_roi*i+j] < SPNumber ) {
                adjacent_sp[size_roi*i+count_adjacent[i]] = tmp_adjacent_sp[size_roi*i+j];
                count_adjacent[i]++;

            }

        }

    }

    delete[] tmp_adjacent_sp;
    delete[] tmp_count_adjacent;

}

void IBIS::mean_seeds() {

    for( int i=0; i<SPNumber; i++ ) {
        inv[ i ] = 1.f / countPx[ i ];
    }

    for( int i=0; i<SPNumber; i++ ) {
        if( countPx[ i ] > 0 ) {
            Xseeds[ i ] = Xseeds_Sum[ i ] * inv[ i ];
            Yseeds[ i ] = Yseeds_Sum[ i ] * inv[ i ];
            lseeds[ i ] = lseeds_Sum[ i ] * inv[ i ];
            aseeds[ i ] = aseeds_Sum[ i ] * inv[ i ];
            bseeds[ i ] = bseeds_Sum[ i ] * inv[ i ];

        }

    }

}

float IBIS::get_complexity() {
    int count_px_processed = 0;
    for( int i=0; i<size; i++ ) {
        if( processed[ i ] )
            count_px_processed++;

    }

    return float(count_px_processed) / float(size);

}

double IBIS::now_ms(void)
{
    double milliseconds_since_epoch = (double) (std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
    return milliseconds_since_epoch;

}

// return the current actual SP number (the method may reduce the actual Sp number).
void IBIS::enforceConnectivity()
{
    //local var
    int label = 0;
    int i, j, k;
    int n, c, count;
    int x, y;
    int ind;
    int oindex, adjlabel;
    int nindex;
    const int dx4[4] = { -1,  0,  1,  0 };
    const int dy4[4] = { 0, -1,  0,  1 };
    int* nlabels = new int[ size ];
    std::fill( nlabels, nlabels + size, -1 );

    oindex = 0;
    adjlabel = 0;

    for (j = 0; j < height; j++)
    {
        for (k = 0; k < width; k++)
        {
            if (nlabels[oindex] < 0)
            {
                nlabels[oindex] = label;// !! labels[oindex] --> label

                x_vec[0] = k;
                y_vec[0] = j;

                for (n = 0; n < 4; n++)
                {
                    x = x_vec[0] + dx4[n];
                    y = y_vec[0] + dy4[n];

                    if ((x >= 0 && x < width) && (y >= 0 && y < height))
                    {
                        nindex = y*width + x;

                        if (nlabels[nindex] >= 0)
                            adjlabel = nlabels[nindex];
                    }
                }

                count = 1;
                for (c = 0; c < count; c++)
                {
                    for (n = 0; n < 4; n++)
                    {
                        x = x_vec[c] + dx4[n];
                        y = y_vec[c] + dy4[n];
                        if ((x >= 0 && x < width) && (y >= 0 && y < height))
                        {
                            nindex = y*width + x;

                            if (nlabels[nindex] < 0 && labels[oindex] == labels[nindex])
                            {
                                x_vec[count] = x;
                                y_vec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }

                if (count <= minSPSizeThreshold)
                {
                    for (c = 0; c < count; c++)
                    {
                        ind = y_vec[c] * width + x_vec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }

    for (i = 0; i < size; i++)
        labels[i] = nlabels[i];

    delete[] nlabels;

}

void IBIS::init() {
    int ii, step;

    // image lab buffer
    avec = new float[size];
    bvec = new float[size];
    lvec = new float[size];

    //store mean distance between 2 seeds info
    SPTypicalLength = (int)(std::sqrt((float)(size) / (float)(maxSPNumber))) + 1;

    // compacity weight
    invwt = (float)SPTypicalLength / compacity;
    invwt = 1.0f / (invwt * invwt);

    // inferior limit for superpixels size
    minSPSizeThreshold = (size / maxSPNumber) / 4;

    // set the top mask size
    index_mask = 1;
    while( SPTypicalLength > ( pow( 2.0, index_mask ) + 1 ) )
        index_mask++;
    index_mask--;

    // set mask size
    mask_size = new int[ index_mask ];
    for( int k=0; k<index_mask; k++ )
        mask_size[ k ] = pow( 2.0, k+1 ) + 1;

    // center of the first block
    start_xy = mask_size[ index_mask - 2 ] - 1;

    // enforce connectivity buffer
    x_vec = new int[ size ];
    y_vec = new int[ size ];

    // visu the processed pixels
    processed = new int[ size ];
    std::fill( processed, processed + size, 0 );

    // precalculate vertical indexes
    vertical_index = new int[ height + mask_size[ index_mask - 1 ] ];
    for( int k=0; k < height + mask_size[ index_mask - 1 ]; k++ ) {
        vertical_index[ k ] = k*width;
    }

    // output labels buffer
    labels = new int[size];

    // repartition of pixels at start
    initial_repartition = new int[size];

    // limit for masks definittion
    y_limit = height + mask_size[ index_mask-1 ];
    x_limit = width + mask_size[ index_mask-1 ];

    // create MASKs
    step = mask_size[ index_mask-1 ];
    ii=0;
    for( int y=start_xy; y<y_limit; y+=step ) {
        for( int x=start_xy; x<x_limit; x+=step ) {
            ii++;

        }

    }

    count_mask = ii;
    mask_buffer = new MASK[ count_mask ];
    ii=0;
    for( int y=start_xy; y<y_limit; y+=step ) {
        for( int x=start_xy; x<x_limit; x+=step ) {
            mask_buffer[ ii ].init( pow(2.0, index_mask+1), y, x, index_mask-1, this );
            ii++;

        }

    }

}

void IBIS::reset() {
    int index_xy;

    st4 = 0;
    st2 = 0;
    st3 = 0;

    std::fill( countPx, countPx + maxSPNumber, 0 );
    std::fill( labels, labels + size, -1 );

    for( int i=0; i < SPNumber; i++ ) {
        Xseeds[ i ] = Xseeds_init[ i ];
        Yseeds[ i ] = Yseeds_init[ i ];

        index_xy = vertical_index[ (int) Yseeds[ i ] ] + Xseeds[ i ];

        lseeds[ i ] = lvec[ index_xy ];
        aseeds[ i ] = avec[ index_xy ];
        bseeds[ i ] = bvec[ index_xy ];
    }

    memset( lseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( aseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( bseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Xseeds_Sum, 0, sizeof( float ) * maxSPNumber );
    memset( Yseeds_Sum, 0, sizeof( float ) * maxSPNumber );

    for( int i=0; i<count_mask; i++ )
        mask_buffer[ i ].reset();

}

void IBIS::getLAB( cv::Mat* img ) {
    cv::Mat lab_image;
    cv::cvtColor(*img, lab_image, cv::COLOR_BGR2Lab, 0);

    int ii = 0;
    for (int i = 0; i < size * 3; i += 3) {
#if MATLAB_lab
        lvec[ii] = lab_image.ptr()[i] * 100 / 255;
        avec[ii] = lab_image.ptr()[i + 1] - 128;
        bvec[ii] = lab_image.ptr()[i + 2] - 128;
#else
        lvec[ii] = lab_image.ptr()[i];
        avec[ii] = lab_image.ptr()[i + 1];
        bvec[ii] = lab_image.ptr()[i + 2];
#endif
        ii++;
    }

}

void IBIS::process( cv::Mat* img ) {
#if OUTPUT_log
    double lap;
    double plab;
#endif

    if( size == 0 ) {
        size = img->cols * img->rows;
        width = img->cols;
        height = img->rows;

        // initialise all the buffer and inner parameters
        init();

        // STEP 1 : initialize with fix grid seeds value
        initSeeds();

    }

    // convert to Lab
#if OUTPUT_log
    plab = now_ms();
#endif

    getLAB( img );

    // prepare value to compute a picture
    reset();

#if OUTPUT_log
    plab = now_ms() - plab;
    lap = now_ms();
#endif

    // STEP 2 : process IBIS
    mask_propagate_SP();

#if OUTPUT_log
    st3 = now_ms() - lap;
#endif

    // STEP 3 : post processing
#if OUTPUT_log
    lap = now_ms();
#endif

    enforceConnectivity();

#if OUTPUT_log
    st4 = now_ms() - lap;
#endif

    // output log
#if OUTPUT_log
    printf("-----------------\n");
    printf("PERF_T %lf\n", st3+st4+plab);
    printf("pre-processing\t\t%lf\t ms\n", plab);
    printf("IBIS.process\t\t%lf\t ms\n", st3);
    printf("IBIS.post_process\t%lf\t ms\n", st4);

    #if MASK_chrono
    float chrono[4] = { 0.f };
    for( int i=0; i < count_mask; i++ )
        mask_buffer[i].get_chrono( chrono );

    float total_chrono = chrono[0] + chrono[1] + chrono[2] + st2;

    printf("-----------------------------------\n");
    printf("Pixels processed:\t%lf\t\%\n", get_complexity()*100 );
    printf("-----------------\n");
    printf("\tMASK.angular_assign()\t\t%lf \%\n", chrono[0]/total_chrono);
    printf("\tMASK.fill_mask()\t\t%lf \%\n", chrono[1]/total_chrono);
    printf("\tMASK.assign_last()\t\t%lf \%\n", chrono[2]/total_chrono);
    printf("\tIBIS.mean_seeds()\t\t%lf \%\n", st2/total_chrono);

    #if THREAD_count > 1
    printf("-----------------\n");
    printf("multi-thread accel:\t\t\t%lf times\n", total_chrono/st3);
    printf("-----------------\n");
    #endif

    #endif

#endif

}

void IBIS::mask_propagate_SP() {

    double lap;
    st2=0;

    for( int mask_index=index_mask-1; mask_index>=0; mask_index--) {

#if THREAD_count > 1
#pragma omp parallel for num_threads(THREAD_count)
#endif
        for( int i=0; i<count_mask; i++ ) {
            mask_buffer[ i ].process( mask_index );

        }

#if MASK_chrono
        lap = now_ms();
#endif
        mean_seeds();
#if MASK_chrono
        st2 += now_ms() - lap;
#endif

#if VISU
        imagesc( std::string("processed"), processed, width, height );
        imagesc( std::string("labels"), labels, width, height );
        cv::waitKey(0);
#endif

    }

}

// ------------------------------------------------------------------------------- IBIS::MASK

void IBIS::MASK::init( int size_in, int y_in, int x_in, int mask_level, IBIS* ptr_IBIS ) {
    IBIS_data = ptr_IBIS;

    x_var = new int[ size_in ];
    y_var = new int[ size_in ];
    xy_var = new int[ size_in ];
    limit_value_fill = new int[4];
    count_last=0;

    size = size_in;
    mask_index = mask_level;

    x = x_in;
    y = y_in;

    filled = false;
    int limit_value = ( IBIS_data->mask_size[ mask_index ] - 1 )/2;
    limit_value_fill[0] = ( x - limit_value < 0 ) ? 0 : x - limit_value;
    limit_value_fill[1] = ( x + limit_value > IBIS_data->width-1 ) ? IBIS_data->width-1 : x + limit_value;
    limit_value_fill[2] = ( y - limit_value < 0 ) ? 0 : y - limit_value;
    limit_value_fill[3] = ( y + limit_value >= IBIS_data->height ) ? IBIS_data->height - 1 : y + limit_value;

    if( mask_index > 0 ) {
        //generate sub_masks
        sub_mask = new IBIS::MASK[ 4 ];

        sub_mask[0].init( pow(2.0, mask_index+1), y - pow(2.0, mask_index-1), x - pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[1].init( pow(2.0, mask_index+1), y - pow(2.0, mask_index-1), x + pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[2].init( pow(2.0, mask_index+1), y + pow(2.0, mask_index-1), x - pow(2.0, mask_index-1), mask_index-1, IBIS_data );
        sub_mask[3].init( pow(2.0, mask_index+1), y + pow(2.0, mask_index-1), x + pow(2.0, mask_index-1), mask_index-1, IBIS_data );

    }
    else {
        last_parent = new int[4];
        last_px_x = new int[9];
        last_px_y = new int[9];
        last_px_xy = new int[9];

        last_parent[0] = -1;
        last_parent[1] = -1;
        last_parent[2] = -1;
        last_parent[3] = -1;

        int j, k, index_y, index_xy;
        for( int index_var_y = -1; index_var_y <= 1; index_var_y++ ) {
            j = y + index_var_y;
            index_y = IBIS_data->vertical_index[ j ];

            for( int index_var_x = -1; index_var_x <= 1; index_var_x++ ) {
                k = x + index_var_x;
                index_xy = index_y + k;

                if( j < IBIS_data->height && k < IBIS_data->width &&
                        !(index_var_x == -1 && index_var_y == -1 ) &&
                        !(index_var_x == -1 && index_var_y == 1 ) &&
                        !(index_var_x == 1 && index_var_y == -1 ) &&
                        !(index_var_x == 1 && index_var_y == 1 ) ) {
                    last_px_x[ count_last ] = k;
                    last_px_y[ count_last ] = j;
                    last_px_xy[ count_last ] = IBIS_data->vertical_index[ j ] + k;

                    count_last++;

                }

            }

        }

    }

    // fill mask coordinates
    generate_mask();

}

IBIS::MASK::~MASK() {
    delete[] x_var;
    delete[] y_var;
    delete[] xy_var;
    delete[] limit_value_fill;

    if( mask_index > 0 ) {
        delete[] sub_mask;

    }
    else {
        delete[] last_px_x;
        delete[] last_px_y;
        delete[] last_px_xy;
        delete[] last_parent;

    }

}

void IBIS::MASK::reset() {
    filled = false;
    Mt1 = 0;
    Mt2 = 0;
    Mt3 = 0;

    if( mask_index > 0 ) {
        sub_mask[ 0 ].reset();
        sub_mask[ 1 ].reset();
        sub_mask[ 2 ].reset();
        sub_mask[ 3 ].reset();

    }
    else {
        last_parent[0] = -1;
        last_parent[1] = -1;
        last_parent[2] = -1;
        last_parent[3] = -1;

    }

}

void IBIS::MASK::get_chrono( float* chrono ) {

    chrono[0]  += Mt1;
    chrono[1]  += Mt2;
    chrono[2]  += Mt3;

    if( mask_index > 0 ) {
        sub_mask[ 0 ].get_chrono( chrono );
        sub_mask[ 1 ].get_chrono( chrono );
        sub_mask[ 2 ].get_chrono( chrono );
        sub_mask[ 3 ].get_chrono( chrono );

    }

}

void IBIS::MASK::assign_labels( int y, int x, int index_xy, int value ) {
    IBIS_data->labels[ index_xy ] = value;

    IBIS_data->countPx[ value ]++;
    IBIS_data->lseeds_Sum[ value ]+=IBIS_data->lvec[ index_xy ];
    IBIS_data->aseeds_Sum[ value ]+=IBIS_data->avec[ index_xy ];
    IBIS_data->bseeds_Sum[ value ]+=IBIS_data->bvec[ index_xy ];
    IBIS_data->Xseeds_Sum[ value ]+=x;
    IBIS_data->Yseeds_Sum[ value ]+=y;

}

void IBIS::MASK::generate_mask() {
    int* tmp_x_var = new int[ size ];
    int* tmp_y_var = new int[ size ];
    int limit_val, vertical_val, value_assign;
    int k = mask_index;
    int table_index = 0;

    if( k > 0 ) {
        limit_val = IBIS_data->mask_size[ k - 1 ] - 1;
        vertical_val = -limit_val;

        for( int index_var=1; index_var <= IBIS_data->mask_size[ k ]; index_var++ ) {
            if( index_var == 1 ) { //top border

                value_assign = -limit_val;
                for( int i=table_index; i<=table_index+limit_val; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + vertical_val;

                    value_assign += 2;
                }

                table_index += limit_val + 1;
                vertical_val += 2;

            }
            else if( index_var > 1 && index_var < IBIS_data->mask_size[ k - 1 ] ) { // vertical border

                value_assign = -limit_val;
                for( int i=table_index; i<table_index+2; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + vertical_val;

                    value_assign = limit_val;
                }

                table_index += 2;
                vertical_val += 2;

            }
            else { // bot border

                value_assign = -limit_val;
                for( int i=table_index; i<=table_index+limit_val; i++ ) {
                    tmp_x_var[ i ] = x + value_assign;
                    tmp_y_var[ i ] = y + limit_val;

                    value_assign += 2;
                }

            }

        }

        // eliminate impossible coordinates
        count_var = 0;
        for( int i=0; i <= table_index + limit_val; i++ ) {
            if( tmp_x_var[ i ] >= 0 && tmp_x_var[ i ] < IBIS_data->width && tmp_y_var[ i ] >= 0 && tmp_y_var[ i ] < IBIS_data->height  ) {
                x_var[ count_var ] = tmp_x_var[ i ];
                y_var[ count_var ] = tmp_y_var[ i ];
                xy_var[ count_var ] = IBIS_data->vertical_index[ tmp_y_var[ i ] ] + tmp_x_var[ i ];
                count_var++;

            }

        }

    }
    else {
        tmp_x_var[ 0 ] = x-1;
        tmp_y_var[ 0 ] = y-1;

        tmp_x_var[ 1 ] = x+1;
        tmp_y_var[ 1 ] = y-1;

        tmp_x_var[ 2 ] = x-1;
        tmp_y_var[ 2 ] = y+1;

        tmp_x_var[ 3 ] = x+1;
        tmp_y_var[ 3 ] = y+1;

        count_var = 0;
        for( int i=0; i < 4; i++ ) {
            if( tmp_x_var[ i ] >= 0 && tmp_x_var[ i ] < IBIS_data->width && tmp_y_var[ i ] >= 0 && tmp_y_var[ i ] < IBIS_data->height  ) {
                x_var[ count_var ] = tmp_x_var[ i ];
                y_var[ count_var ] = tmp_y_var[ i ];
                xy_var[ count_var ] = IBIS_data->vertical_index[ tmp_y_var[ i ] ] + tmp_x_var[ i ];
                count_var++;

            }

        }

    }

    if( count_var == 0 )
        filled = 1;

    delete[] tmp_x_var;
    delete[] tmp_y_var;

}

void IBIS::MASK::process( int mask_level) {
    bool stat;
    double lap;

    if( ! filled ) {
        if( mask_level == mask_index ) {
#if MASK_chrono
            lap = IBIS_data->now_ms();
#endif
            stat = angular_assign();
#if MASK_chrono
            Mt1 += IBIS_data->now_ms() - lap;
#endif

            if( stat ) {
#if MASK_chrono
                lap = IBIS_data->now_ms();
#endif
                fill_mask();
#if MASK_chrono
                Mt2 += IBIS_data->now_ms() - lap;
#endif

                return;
            }
            else {
                if( mask_index == 0 ) {
#if MASK_chrono
                    lap = IBIS_data->now_ms();
#endif
                    assign_last();
#if MASK_chrono
                    Mt3 += IBIS_data->now_ms() - lap;
#endif

                }

            }

        }
        else {
            if( mask_index > 0 ) {
                sub_mask[ 0 ].process( mask_level );
                sub_mask[ 1 ].process( mask_level );
                sub_mask[ 2 ].process( mask_level );
                sub_mask[ 3 ].process( mask_level );

            }

        }

    }

#if VISU_all
    imagesc( std::string("processed"), IBIS_data->processed, IBIS_data->width, IBIS_data->height );
    cv::waitKey(1);
#endif

}

int IBIS::MASK::assign_last_px( int y, int x, int index_xy ) {
    int best_sp = -1;
    float dist_xy;
    float l, a, b;
    float dist_lab;
    int index_sp;
    float D=-1.f;
    float total_dist;

    for( int i=0; i<count_last_parent; i++ ) {
        index_sp = last_parent[ i ];

        l = IBIS_data->lvec[ index_xy ];
        a = IBIS_data->avec[ index_xy ];
        b = IBIS_data->bvec[ index_xy ];

        dist_lab = ( l - IBIS_data->lseeds[ index_sp ]) * ( l - IBIS_data->lseeds[ index_sp ]) +
                   ( a - IBIS_data->aseeds[ index_sp ]) * ( a - IBIS_data->aseeds[ index_sp ]) +
                   ( b - IBIS_data->bseeds[ index_sp ]) * ( b - IBIS_data->bseeds[ index_sp ]);

        dist_xy = ( x - IBIS_data->Xseeds[ index_sp ] ) * ( x - IBIS_data->Xseeds[ index_sp ] ) +
                  ( y - IBIS_data->Yseeds[ index_sp ] ) * ( y - IBIS_data->Yseeds[ index_sp ] );


        total_dist = dist_lab + dist_xy * IBIS_data->invwt;

        if( total_dist < D || D < 0) {
            best_sp = index_sp;
            D = total_dist;

        }

    }

    return best_sp;

}

int IBIS::MASK::assign_px( int y, int x, int index_xy ) {
    int best_sp = -1;
    float dist_xy;
    float l, a, b;
    float dist_lab;
    int index_sp;
    float D=-1.f;
    float total_dist;

    for( int i=0; i<IBIS_data->count_adjacent[IBIS_data->initial_repartition[index_xy]]; i++ ) {
        index_sp = IBIS_data->adjacent_sp[ size_roi*IBIS_data->initial_repartition[index_xy] + i ];

        l = IBIS_data->lvec[ index_xy ];
        a = IBIS_data->avec[ index_xy ];
        b = IBIS_data->bvec[ index_xy ];

        dist_lab = ( l - IBIS_data->lseeds[ index_sp ]) * ( l - IBIS_data->lseeds[ index_sp ]) +
                   ( a - IBIS_data->aseeds[ index_sp ]) * ( a - IBIS_data->aseeds[ index_sp ]) +
                   ( b - IBIS_data->bseeds[ index_sp ]) * ( b - IBIS_data->bseeds[ index_sp ]);

        dist_xy = ( x - IBIS_data->Xseeds[ index_sp ] ) * ( x - IBIS_data->Xseeds[ index_sp ] ) +
                  ( y - IBIS_data->Yseeds[ index_sp ] ) * ( y - IBIS_data->Yseeds[ index_sp ] );


        total_dist = dist_lab + dist_xy * IBIS_data->invwt;

        if( total_dist < D || D < 0) {
            best_sp = index_sp;
            D = total_dist;

        }

    }

    return best_sp;
}

void IBIS::MASK::assign_last() {

    int value;

    for( int index_last=0; index_last < count_last; index_last++ ) {
        if( IBIS_data->labels[ last_px_xy[ index_last ] ] < 0 ) {
            value = assign_last_px( last_px_y[ index_last ], last_px_x[ index_last ], last_px_xy[ index_last ] );
            assign_labels( last_px_y[ index_last ], last_px_x[ index_last ], last_px_xy[ index_last ], value );

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ last_px_xy[ index_last ] ] = 1;
#endif

        }

    }

}

void IBIS::MASK::fill_mask() {

    int j;
    int index_y;

    for( int index_var_y = limit_value_fill[2]; index_var_y <= limit_value_fill[3]; index_var_y++ ) {
        for( int index_var_x = limit_value_fill[0]; index_var_x <= limit_value_fill[1]; index_var_x++ ) {
            assign_labels( index_var_y, index_var_x, IBIS_data->vertical_index[ index_var_y ] + index_var_x, angular );

        }
        
    }

    filled = true;

}

bool IBIS::MASK::angular_assign() {
    int j, k;
    int best_sp;
    int index_xy;
    bool first=true;
    bool output = true;
    int ref_sp;
    angular = -1;

    if( mask_index > 0 ) {
        for( int index_var=0; index_var < count_var; index_var++ ) {
            j = y_var[ index_var ];
            k = x_var[ index_var ];
            index_xy = xy_var[ index_var ];

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ index_xy ] = 1;
#endif

            if( IBIS_data->labels[ index_xy ] < 0 ) {
                // assign px
                angular = assign_px( j, k, index_xy );

            }
            else
                angular = IBIS_data->labels[ index_xy ];

            if( first ) {
                first = false;
                ref_sp = angular;

            }

            if( ref_sp != angular && !first ) {
                 return false;

            }

        }

    }
    else {
        count_last_parent = 0;

        for( int index_var=0; index_var < count_var; index_var++ ) {
            j = y_var[ index_var ];
            k = x_var[ index_var ];
            index_xy = xy_var[ index_var ];

#if VISU || VISU_all || MASK_chrono
            IBIS_data->processed[ index_xy ] = 1;
#endif

            if( IBIS_data->labels[ index_xy ] < 0 ) {
                // assign px
                best_sp = assign_px( j, k, index_xy );

                //IBIS_data->labels[ index_xy ] = best_sp;
                assign_labels( j, k, IBIS_data->vertical_index[ j ] + k, best_sp );

            }

            angular = IBIS_data->labels[ index_xy ];

            if( (last_parent[0] != angular) &&
                (last_parent[1] != angular) &&
                (last_parent[2] != angular) &&
                (last_parent[3] != angular) ) {

                last_parent[ count_last_parent ] = angular;
                count_last_parent++;

            }

            if( first ) {
                ref_sp = angular;
                first = false;

            }

            if( ref_sp != angular && !first )
                output = false;

        }

    }

    return output;

}
