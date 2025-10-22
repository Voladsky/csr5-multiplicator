#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct {
	uint32_t* 	rows;
	uint32_t* 	cols;
	double*	data;
	uint32_t 	nnz;
	uint32_t 	nrows;
	uint32_t 	ncols;
} COOMat;

typedef struct {
	uint32_t* 	row_ptr;
	uint32_t* 	col_idx;
	double*	data;
	uint32_t 	nnz;
	uint32_t 	nrows;
	uint32_t 	ncols;
} CSRMat;

typedef struct {
	char* bit_flag;
	uint32_t* y_offset;
	uint32_t* seg_offset;
	uint32_t* empty_offset;
} TileDesc;


typedef struct {
	uint32_t* 		row_ptr;
	uint32_t*		col_idx;
	double*		data;
	uint32_t* 	tile_ptr;
	TileDesc* 	tile_desc;
	uint32_t 		nnz;
	uint32_t 		nrows;
	uint32_t 		ncols;
	uint32_t 		omega;
	uint32_t 		sigma;
	uint32_t 		p;
} CSR5Mat;

CSRMat coo2csr(COOMat coo) {
	CSRMat csr;
	csr.nnz 	= coo.nnz;
	csr.nrows 	= coo.nrows;
	csr.ncols 	= coo.ncols;
	csr.data 	= malloc(csr.nnz * sizeof(double));
	csr.col_idx = malloc(csr.nnz * sizeof(uint32_t));
	csr.row_ptr = calloc((csr.nrows + 1), sizeof(*csr.row_ptr));
	
	memcpy(csr.data, coo.data, coo.nnz * sizeof(double));
	memcpy(csr.col_idx, coo.cols, coo.nnz * sizeof(uint32_t));
	
	for (uint32_t i = 0; i < csr.nnz; ++i) {
		csr.row_ptr[coo.rows[i] + 1]++;
	}
	for (uint32_t i = 1; i <= csr.nrows; ++i) {
		csr.row_ptr[i] += csr.row_ptr[i - 1];
	}
	
	return csr;
}

uint32_t cmp(const void* arg1, const void* arg2) {
	uint32_t a = *(uint32_t*)arg1;
	uint32_t b = *(uint32_t*)arg2;
	return (a > b) - (a < b);
}

uint32_t remove_sign(uint32_t x) {
	return x &= 0x7FFFFFFF;
}

void segsum(uint32_t* in, uint32_t* flag, uint32_t length) {
	for (uint32_t i = 0; i < length; ++i) {
		if (flag[i]) {
			uint32_t j = i + 1;
			while (!flag[j] && j < length) {
				in[i] += in[j];
				++j;
			}
		}
		else {
			in[i] = 0;
		}
	}
}

void exclusive_scan(uint32_t* in, uint32_t length) {
	uint32_t tmp = in[0];
	in[0] = 0;
	for (uint32_t i = 1; i < length; ++i) {
		uint32_t cur = in[i];
		in[i] = in[i - 1] + tmp;
		tmp = cur;
	}
}

uint32_t reduction_sum(uint32_t* in, uint32_t length) {
	uint32_t sum = 0;
	for (uint32_t i = 0; i < length; ++i) {
		sum += in[i];
	}
	return sum;
}

uint32_t binary_search(const uint32_t* arr, uint32_t size, uint32_t value) {
	uint32_t left = 0, right = size - 1;
	while (left <= right) {
		uint32_t mid = left + (right - left) / 2;
		if (arr[mid] <= value) {
			left = mid + 1;
		}
		else {
			right = mid - 1;
		}
	}
	return left - 1; // returns the largest index where arr[index] <= value
}

bool is_row_start(const uint32_t* row_ptr, uint32_t nrows, uint32_t pos) {
	for (uint32_t r = 0; r < nrows; ++r) {
		if (row_ptr[r] == pos) {
			return true;
		}
	}
	return false;
}

CSR5Mat csr2csr5(CSRMat csr, uint32_t omega, uint32_t sigma) {
	uint32_t p = (uint32_t)ceil(csr.nnz * 1.0 / (sigma * omega));
	CSR5Mat csr5;
	csr5.nnz		= csr.nnz;
	csr5.nrows		= csr.nrows;
	csr5.ncols		= csr.ncols;
	csr5.omega		= omega;
	csr5.sigma 		= sigma;
	csr5.p			= p;
	csr5.data		= malloc(csr5.nnz * sizeof(*csr5.data));
	csr5.col_idx 	= malloc(csr5.nnz * sizeof(*csr5.col_idx));
	csr5.row_ptr 	= malloc((csr5.nrows + 1) * sizeof(*csr5.row_ptr));
	csr5.tile_ptr 	= malloc((p + 1) * sizeof(*csr5.tile_ptr));
	csr5.tile_desc	= malloc((p + 1) * sizeof(*csr5.tile_desc));
	
	memcpy(csr5.row_ptr, csr.row_ptr, (csr5.nrows + 1) * sizeof(uint32_t));
	memcpy(csr5.col_idx, csr.col_idx, csr.nnz * sizeof(uint32_t));
    memcpy(csr5.data, csr.data, csr.nnz * sizeof(double));
	
	//csr5.bit_flag 		= malloc(omega * sigma * sizeof(uint32_t));
	//csr5.y_offset 		= malloc(omega * sizeof(uint32_t));
	//csr5.seg_offset 	= malloc(omega * sizeof(uint32_t));
	// empty offset is allocated later

	// generate tile_ptr
	for (uint32_t tid = 0; tid < p + 1; ++tid) {
		uint32_t bnd = tid * omega * sigma;
		csr5.tile_ptr[tid] = binary_search(csr5.row_ptr, csr5.nrows + 1, bnd);
	}
	
	for (uint32_t tid = 0; tid < p; ++tid) {
		for (uint32_t rid = csr5.tile_ptr[tid]; rid <= csr5.tile_ptr[tid + 1]; ++rid) {
			if (csr5.row_ptr[rid] == csr5.row_ptr[rid + 1]) {
				csr5.tile_ptr[tid] |= 0x80000000; // inverse first bit
				break;
			}
		}
	}
	
	for (uint32_t tid = 0; tid < p; ++tid) {
		
		TileDesc* desc		= &csr5.tile_desc[tid];
        desc->bit_flag 		= calloc(omega * sigma, sizeof(char));
        desc->y_offset 		= malloc(omega * sizeof(uint32_t));
        desc->seg_offset 	= malloc(omega * sizeof(uint32_t));
		

		uint32_t tile_start = tid * omega * sigma;
		uint32_t tile_end	= (tid == p - 1) ? csr5.nnz : tile_start + omega * sigma;
		uint32_t tile_size	= tile_end - tile_start;

		// generate bit_flag
		if (tile_size > 0) {
			desc->bit_flag[0] = 1;
		}

		// Set bit_flag for row starts within this tile
		for (uint32_t pos = tile_start; pos < tile_end; ++pos) {
			if (is_row_start(csr5.row_ptr, csr5.nrows, pos)) {
				uint32_t tile_local_pos = pos - tile_start;
				if (tile_local_pos < omega * sigma) {
					desc->bit_flag[tile_local_pos] = 1;
				}
			}
		}
		
		// generate y_offset and seg_offset
		uint32_t* tmp_bit = calloc(omega, sizeof(uint32_t));
		for (uint32_t i = 0; i < omega; ++i) {
			desc->y_offset[i] = 0;
			for (uint32_t j = 0; j < sigma; ++j) {
				desc->y_offset[i] += desc->bit_flag[i * omega + j];
				tmp_bit[i] 	|= desc->bit_flag[i * omega + j];
			}
			desc->seg_offset[i] = 1 - tmp_bit[i];
		}
		exclusive_scan(desc->y_offset, omega);
		segsum(desc->seg_offset, tmp_bit, omega);
		
		// generate empty_offset
		if (csr5.tile_ptr[tid] >> 31) {
			uint32_t length 			= reduction_sum(desc->bit_flag, omega * sigma);
			desc->empty_offset 	= malloc(length * sizeof(uint32_t));
			uint32_t eid = 0;
			for (uint32_t i = 0; i < omega; ++i) {
				for (uint32_t j = 0; j < sigma; ++j) {
					if (desc->bit_flag[i * omega + j]) {
						uint32_t ptr = tid * omega * sigma + i * sigma + j;
						uint32_t idx = binary_search( csr5.row_ptr, csr5.nrows + 1, ptr) - 1;
						idx = idx - remove_sign(csr5.tile_ptr[tid]);
						desc->empty_offset[eid] = idx + 1;
						++eid;
					}
				}
			}
		}
		else {
			desc->empty_offset = NULL;
		}
		free(tmp_bit);
	}
	
	return csr5;
}

void printcoo(COOMat coo) {
	for (uint32_t i = 0; i < coo.nnz; ++i) {
		printf("%d, %d, %.2f\n", coo.rows[i], coo.cols[i], coo.data[i]);
	}
}
void printcsr(CSRMat csr) {
	printf("data:\t\t");
	for (uint32_t i = 0; i < csr.nnz; ++i) {
		printf("%.2f ", csr.data[i]);
	}
	printf("\ncol_idx:\t");
	for (uint32_t i = 0; i < csr.nnz; ++i) {
		printf("%d ", csr.col_idx[i]);
	}
	printf("\nrow_ptr:\t");
	for (uint32_t i = 0; i < csr.nrows + 1; ++i) {
		printf("%d ", csr.row_ptr[i]);
	}
	printf("\n");
}

void printcsr5(CSR5Mat csr5) {
    printf("=== CSR5 Matrix Format ===\n");
    printf("Dimensions: %d x %d\n", csr5.nrows, csr5.ncols);
    printf("Nonzeros: %d\n", csr5.nnz);
    printf("Parameters: omega=%d, sigma=%d, p=%d\n", csr5.omega, csr5.sigma, csr5.p);
    printf("\n");
    
    // Print basic CSR arrays
    printf("1. Basic CSR Arrays:\n");
    printf("   row_ptr (%d elements): ", csr5.nrows + 1);
    for (uint32_t i = 0; i <= csr5.nrows && i < 20; ++i) {
        printf("%d ", csr5.row_ptr[i]);
        if (i == 10 && csr5.nrows > 20) {
            printf("... ");
            i = csr5.nrows - 5;
        }
    }
    printf("\n");
    
    printf("   col_idx (%d elements): ", csr5.nnz);
    for (uint32_t i = 0; i < csr5.nnz && i < 20; ++i) {
        printf("%d ", csr5.col_idx[i]);
        if (i == 10 && csr5.nnz > 20) {
            printf("... ");
            i = csr5.nnz - 5;
        }
    }
    printf("\n");
    
    printf("   data (%d elements): ", csr5.nnz);
    for (uint32_t i = 0; i < csr5.nnz && i < 10; ++i) {
        printf("%.1f ", csr5.data[i]);
        if (i == 5 && csr5.nnz > 10) {
            printf("... ");
            i = csr5.nnz - 3;
        }
    }
    printf("\n\n");
    
    // Print tile pointers
    printf("2. Tile Pointers (p+1=%d elements):\n", csr5.p + 1);
    for (uint32_t tid = 0; tid <= csr5.p; ++tid) {
        if (csr5.tile_ptr[tid] & 0x80000000) {
            printf("   Tile %d: %d (has empty rows)\n", tid, remove_sign(csr5.tile_ptr[tid]));
        } else {
            printf("   Tile %d: %d\n", tid, csr5.tile_ptr[tid]);
        }
    }
    printf("\n");
    
    // Print tile descriptors for first few tiles
    uint32_t tiles_to_print = (csr5.p < 3) ? csr5.p : 3;
    printf("3. Tile Descriptors (showing first %d tiles):\n", tiles_to_print);
    
    for (uint32_t tid = 0; tid < tiles_to_print; ++tid) {
        TileDesc* desc = &csr5.tile_desc[tid];
        printf("   Tile %d:\n", tid);
        
        // Print bit_flag
        printf("     bit_flag (%dx%d): ", csr5.omega, csr5.sigma);
        uint32_t elements_to_print = (csr5.omega * csr5.sigma < 16) ? csr5.omega * csr5.sigma : 16;
        for (uint32_t i = 0; i < elements_to_print; ++i) {
            printf("%d", desc->bit_flag[i]);
            if ((i + 1) % csr5.sigma == 0) printf(" "); // New line every sigma elements
        }
        if (elements_to_print < csr5.omega * csr5.sigma) printf("...");
        printf("\n");
        
        // Print y_offset
        printf("     y_offset (%d): [", csr5.omega);
        for (uint32_t i = 0; i < csr5.omega; ++i) {
            printf("%d", desc->y_offset[i]);
            if (i < csr5.omega - 1) printf(", ");
        }
        printf("]\n");
        
        // Print seg_offset
        printf("     seg_offset (%d): [", csr5.omega);
        for (uint32_t i = 0; i < csr5.omega; ++i) {
            printf("%d", desc->seg_offset[i]);
            if (i < csr5.omega - 1) printf(", ");
        }
        printf("]\n");
        
        // Print empty_offset if exists
        if (desc->empty_offset != NULL) {
            uint32_t empty_len = 0;
            // Calculate length by counting TRUEs in bit_flag
            for (uint32_t i = 0; i < csr5.omega * csr5.sigma; ++i) {
                if (desc->bit_flag[i]) empty_len++;
            }
            printf("     empty_offset (%d): [", empty_len);
            for (uint32_t i = 0; i < empty_len && i < 10; ++i) {
                printf("%d", desc->empty_offset[i]);
                if (i < empty_len - 1 && i < 9) printf(", ");
            }
            if (empty_len > 10) printf("...");
            printf("]\n");
        } else {
            printf("     empty_offset: NULL\n");
        }
        printf("\n");
    }
    
    if (csr5.p > tiles_to_print) {
        printf("   ... and %d more tiles\n", csr5.p - tiles_to_print);
    }
    
    // Print memory usage summary
    printf("4. Memory Usage Summary:\n");
    size_t basic_size = (csr5.nrows + 1) * sizeof(uint32_t) + 
                       csr5.nnz * sizeof(uint32_t) + 
                       csr5.nnz * sizeof(double);
    
    size_t tile_ptr_size = (csr5.p + 1) * sizeof(uint32_t);
    
    size_t tile_desc_size = 0;
    for (uint32_t tid = 0; tid < csr5.p; ++tid) {
        TileDesc* desc = &csr5.tile_desc[tid];
        tile_desc_size += csr5.omega * csr5.sigma * sizeof(char) +  // bit_flag
                         csr5.omega * sizeof(uint32_t) * 2;              // y_offset + seg_offset
        
        if (desc->empty_offset != NULL) {
            uint32_t empty_len = 0;
            for (uint32_t i = 0; i < csr5.omega * csr5.sigma; ++i) {
                if (desc->bit_flag[i]) empty_len++;
            }
            tile_desc_size += empty_len * sizeof(uint32_t);
        }
    }
    
    size_t total_size = basic_size + tile_ptr_size + tile_desc_size;
    
    printf("   Basic CSR arrays: %zu bytes\n", basic_size);
    printf("   Tile pointers: %zu bytes\n", tile_ptr_size);
    printf("   Tile descriptors: %zu bytes\n", tile_desc_size);
    printf("   TOTAL: %zu bytes (%.2f KB)\n", total_size, total_size / 1024.0);
    printf("   Overhead: %.1f%%\n", 
           (tile_ptr_size + tile_desc_size) * 100.0 / basic_size);
}

uint32_t main() {
	// Global consts
	uint32_t m = 8, n = 8, nnz = 34;
	uint32_t omega = 4, sigma = 4;
	
	// Define matrix with triplets
	// 1 0 2 3 0 0 4 5
	// 0 1 0 2 0 0 0 0
	// 0 0 0 0 0 0 0 0
	// 1 2 3 4 5 0 6 7
	// 0 1 0 2 0 3 0 0
	// 1 2 0 0 0 0 0 0
	// 0 1 2 3 4 5 6 7
	// 1 2 3 4 5 6 7 8
	COOMat coo;
	coo.nnz 	= nnz;
	coo.nrows 	= m;
	coo.ncols 	= n;
	coo.rows	= (uint32_t[])	{ 0, 0, 0, 0, 0, 1, 1, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7 }; 
	coo.cols	= (uint32_t[])	{ 0, 2, 3, 6, 7, 1, 3, 0, 1, 2, 3, 4, 6, 7, 1, 3, 5, 0, 1, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7 };
	coo.data	= (double[]){ 1, 2, 3, 4, 5, 1, 2, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 1, 2, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8 };
	
	CSRMat csr = coo2csr(coo);

	printf("----------------------------\n");
	
	printcoo(coo);
	
	printf("----------------------------\n");
	
	printcsr(csr);
	
	CSR5Mat csr5 = csr2csr5(csr, omega, sigma);
	
	printf("----------------------------\n");
	
	printcsr5(csr5);
	
	
	coo.nnz 	= 7;
	coo.nrows 	= 4;
	coo.ncols 	= 4;
	coo.rows 	= (uint32_t[])	{ 0, 0, 2, 2, 2, 3, 3 };
	coo.cols	= (uint32_t[])	{ 0, 2, 0, 2, 3, 1, 3 };
	coo.data	= (double[]){ 1, 2, 1, 2, 3, 1, 2 };
	
	CSRMat csr2 = coo2csr(coo);
	
	printf("----------------------------\n");
	
	printcsr(csr2);
	
	// Free the slaves of malloc
	free(csr.row_ptr);
	free(csr.col_idx);
	free(csr.data);
	free(csr2.row_ptr);
	free(csr2.col_idx);
	free(csr2.data);
	return 0;
}