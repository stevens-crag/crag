/*
 * SignMatrix.h
 *
 *  Created on: 21.03.2011
 *      Author: Juern
 */

#ifndef SIGNMATRIX_H_
#define SIGNMATRIX_H_

namespace PC
{

/**
 * Data type that stores a matrix containing entries from {-1,0,+1}
 * in an efficient way. In addition the matrix contains a row of
 * column headers and a column of row headers whose types are
 * specified via templates. It can be expanded by inserting rows
 * and columns at any position.
 * It is optimized to make operations on columns somewhat faster
 * than those on rows.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE> class SignMatrix;
struct Row {unsigned int rowid; static const unsigned int UNDEFINED = ~0; Row():rowid(UNDEFINED){}};
struct Col {unsigned int colid; static const unsigned int UNDEFINED = ~0; Col():colid(UNDEFINED){}};
const Row UNDEFINED_ROW;
const Col UNDEFINED_COL;

template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
class SignMatrix
{
private:
	unsigned int		numRows, numCols;
	unsigned int		maxNumRows, maxNumCols;
	int 				firstUnusedRow;
	int 				firstUnusedCol;

	INTERNAL_TYPE*		matrix;
	struct intRowHdr {unsigned int rowno, row; RowHdr rh; int nextUnused; bool used;};
	struct intColHdr {unsigned int colno, col; ColHdr ch; int nextUnused; bool used;};
	std::vector<intRowHdr>	rowHdrs;
	std::vector<intColHdr>	colHdrs;

	static const double	EXTEND_CONDITION = 1.2;
	static const double	EXTEND_FACTOR = 2.0;
	const unsigned int	INTERNAL_SIZE;
	INTERNAL_TYPE		mask1, mask2;

	void				extend(unsigned int rows, unsigned int cols);

	Row					insertRow(unsigned int rowNo, RowHdr rowHdr);
	Col					insertCol(unsigned int colNo, ColHdr colHdr);
	void				free();
public:
						SignMatrix(unsigned int rows = 0, unsigned int cols = 0, unsigned int maxR = 16, unsigned int maxC = 16);
	virtual				~SignMatrix();
	// Access matrix and headers
	unsigned int		getNumRows() const;
	unsigned int		getNumCols() const;
	Row					getRow(unsigned int rowNo) const;
	unsigned int		getRowNumber(Row row) const;
	Col					getCol(unsigned int colNo) const;
	unsigned int		getColNumber(Col col) const;
	Sign				operator()(unsigned int rowNo, unsigned int colNo) const;
	Sign				operator()(Row row, Col col) const;
	Sign				operator()(Row row, unsigned int colNo) const;
	Sign				operator()(unsigned int rowNo, Col col) const;
	void				setEntry(unsigned int rowNo, unsigned int colNo, Sign val);
	void				setEntry(Row row, Col col, Sign val);
	void				setEntry(Row row, unsigned int colNo, Sign val);
	void				setEntry(unsigned int rowNo, Col col, Sign val);
	RowHdr&				rowHdr(unsigned int rowNo);
	RowHdr&				rowHdr(Row col);
	ColHdr&				colHdr(unsigned int colNo);
	ColHdr&				colHdr(Col col);
	// Operations on columns
	bool				columnsEqual(Col col1, Col col2);
	Col					insertZeroCol(unsigned int colNo, ColHdr colHdr = ColHdr());
	Col					insertUnitCol(unsigned int colNo, unsigned int i, ColHdr colHdr = ColHdr());
	Col					insertInvCol(Col sourceCol, unsigned int targetColNo, ColHdr targetColHdr = ColHdr());
	Col					insertColSum(Col sourceCol1, Col sourceCol2, unsigned int targetColNo, ColHdr targetColHdr = ColHdr());
	Col					insertColCopy(Col sourceCol, unsigned int targetColNo);
	Col					insertIntersectionCol(Col sourceCol1, Col sourceCol2, unsigned int targetColNo, ColHdr targetColHdr = ColHdr());
	void				addToCol(Col sourceCol, Col targetCol);
	void				deleteCol(Col col);
	void				moveCol(Col col, unsigned int newColNo);
	// Operations on rows
	Row					insertZeroRow(unsigned int rowNo, RowHdr rowHdr = RowHdr());
	Row					insertRowCopy(Row sourceRow, unsigned int targetRowNo);
	void				deleteRow(Row row);
	void				moveRow(Row row, unsigned int newRowNo);
	// Other matrix operations
	void				clone(const SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>& mat, std::list<Col>& colList);
	void				insertMatrix(const SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>& mat, unsigned int rowOffset, unsigned int colOffset);
	void				printMatrix(std::ostream& os = std::cout) const;
#ifdef DEBUG
	bool				checkInternalIntegrity();
#endif
};


/**
 * Constructor
 *
 * Creates an instance of SignMatrix having the indicated number of
 * rows and columns (initially filled with zeros). If there is an
 * estimate on the final size of the matrix, it can be given in
 * order to avoid later memory reallocations.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::SignMatrix(unsigned int rows, unsigned int cols, unsigned int maxR, unsigned int maxC)
:	numRows(rows),
 	numCols(cols),
	firstUnusedRow(-1),
	firstUnusedCol(-1),
 	matrix(NULL),
 	rowHdrs(),
 	colHdrs(),
	INTERNAL_SIZE(sizeof(INTERNAL_TYPE) << 2)
{
 	maxNumRows = (((std::max(maxR, rows) + INTERNAL_SIZE - 1) / INTERNAL_SIZE) * INTERNAL_SIZE);
 	maxNumCols = std::max(maxC, cols);
	// Initialize masks
	for(unsigned int j = 0; j < INTERNAL_SIZE; j++)
		SignMatrix::mask1 = (SignMatrix::mask1 << 2) | 1/*0b01*/;
	mask2 = mask1 << 1;
	// Allocate memory
	if(maxNumRows * maxNumCols > 0)
	{
		matrix = new INTERNAL_TYPE[maxNumCols * (maxNumRows / INTERNAL_SIZE)];
		assert(matrix != NULL);
		memset(matrix, 0, maxNumCols * (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
	}
	if(maxNumRows > 0)
	{
		rowHdrs.resize(maxNumRows);

		for(unsigned int i = 0; i < numRows; i++)
		{
			rowHdrs[i].used = true;
			rowHdrs[i].row = i;
			rowHdrs[i].rowno = i;
		}
	 	if(numRows < maxNumRows)
	 		firstUnusedRow = numRows;
	 	else
	 		firstUnusedRow = -1;

		for(unsigned int i = numRows; i < maxNumRows; i++)
		{
			rowHdrs[i].used = false;
			rowHdrs[i].nextUnused = i + 1;
		}
		rowHdrs[maxNumRows-1].nextUnused = -1;
	}
	if(maxNumCols > 0)
	{
		colHdrs.resize(maxNumCols);

		for(unsigned int i = 0; i < numCols; i++)
		{
			colHdrs[i].used = true;
			colHdrs[i].col = i;
			colHdrs[i].colno = i;
		}
	 	if(numCols < maxNumCols)
	 		firstUnusedCol = numCols;
	 	else
	 		firstUnusedCol = -1;
		for(unsigned int i = numCols; i < maxNumCols; i++)
		{
			colHdrs[i].used = false;
			colHdrs[i].nextUnused = i + 1;
		}
		colHdrs[maxNumCols-1].nextUnused = -1;
	}
}

/**
 * Destructor
 *
 * Destroys the instance and deallocates any used memory.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::~SignMatrix()
{
	free();
}

/**
 * Cleans up and deallocates any memory used by this object
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::free()
{
	if(matrix != NULL)
	{
		delete [] matrix;
		matrix = NULL;
	}
	rowHdrs.empty();
	colHdrs.empty();

	numRows = 0; numCols = 0;
	maxNumRows = 0; maxNumCols = 0;
}

/**
 * Methods for getting the size of the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline unsigned int SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getNumRows() const
{
	return numRows;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline unsigned int SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getNumCols() const
{
	return numCols;
}

/**
 * Methods for getting the index of a row or column or
 * getting the row or column at a specified index.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline Row SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getRow(unsigned int rowNo) const
{
	assert(rowNo < numRows && rowHdrs[rowHdrs[rowNo].row].used);
	Row row; row.rowid = rowHdrs[rowNo].row;
	return row;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline unsigned int SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getRowNumber(Row row) const
{
	assert(row.rowid < maxNumRows && rowHdrs[row.rowid].rowno < numRows && rowHdrs[row.rowid].used);
	return rowHdrs[row.rowid].rowno;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getCol(unsigned int colNo) const
{
	assert(colNo < numCols && colHdrs[colHdrs[colNo].col].used);
	Col col; col.colid = colHdrs[colNo].col;
	return col;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline unsigned int SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::getColNumber(Col col) const
{
	assert(col.colid < maxNumCols && colHdrs[col.colid].colno < numCols && colHdrs[col.colid].used);
	return colHdrs[col.colid].colno;
}

/**
 * Operators giving access to the entries of the matrix.
 * Index bounds are checked in debug mode only.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline Sign SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::operator()(Row row, Col col) const
{
	assert(row.rowid < maxNumRows && col.colid < maxNumCols && rowHdrs[row.rowid].rowno < numRows && colHdrs[col.colid].colno < numCols && rowHdrs[row.rowid].used && colHdrs[col.colid].used);
	unsigned long long  t = row.rowid + col.colid * maxNumRows;
	return (matrix[t / INTERNAL_SIZE] >> ((t % INTERNAL_SIZE) * 2)) & 3/*0b11*/;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline Sign SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::operator()(unsigned int rowNo, unsigned int colNo) const
{
	assert(rowNo < numRows && colNo < numCols && rowHdrs[rowHdrs[rowNo].row].used && colHdrs[colHdrs[colNo].col].used);
	unsigned long long t = rowHdrs[rowNo].row + colHdrs[colNo].col * maxNumRows;
	return (matrix[t / INTERNAL_SIZE] >> ((t % INTERNAL_SIZE) * 2)) & 3/*0b11*/;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Sign SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::operator()(Row row, unsigned int colNo) const
{
	return operator()(row, getCol(colNo));
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Sign SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::operator()(unsigned int rowNo, Col col) const
{
	return operator()(getRow(rowNo), col);
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::setEntry(Row row, Col col, Sign val)
{
	assert(row.rowid < maxNumRows && col.colid < maxNumCols && rowHdrs[row.rowid].rowno < numRows && colHdrs[col.colid].colno < numCols && rowHdrs[row.rowid].used && colHdrs[col.colid].used);
	unsigned long long t = row.rowid + col.colid * maxNumRows;
	matrix[t / INTERNAL_SIZE] &= ~(3/*0b11*/ << ((t % INTERNAL_SIZE) * 2));
	matrix[t / INTERNAL_SIZE] |= (val & 3/*0b11*/) << ((t % INTERNAL_SIZE) * 2);
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::setEntry(unsigned int rowNo, unsigned int colNo, Sign val)
{
	assert(rowNo < numRows && colNo < numCols && rowHdrs[rowHdrs[rowNo].row].used && colHdrs[colHdrs[colNo].col].used);
	unsigned long long t = rowHdrs[rowNo].row + colHdrs[colNo].col * maxNumRows;
	matrix[t / INTERNAL_SIZE] &= ~(3/*0b11*/ << ((t % INTERNAL_SIZE) * 2));
	matrix[t / INTERNAL_SIZE] |= (val & 3/*0b11*/) << ((t % INTERNAL_SIZE) * 2);
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::setEntry(Row row, unsigned int colNo, Sign val)
{
	setEntry(row, getCol(colNo), val);
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::setEntry(unsigned int rowNo, Col col, Sign val)
{
	setEntry(getRow(rowNo), col, val);
}

/**
 * Operators giving access to the headers of the matrix.
 * Index bounds are checked in debug mode only.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline RowHdr& SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::rowHdr(unsigned int rowNo)
{
	assert(rowNo < numRows && rowHdrs[rowHdrs[rowNo].row].used);
	return rowHdr(getRow(rowNo));
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline RowHdr& SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::rowHdr(Row row)
{
	assert(row.rowid < maxNumRows && rowHdrs[row.rowid].rowno < numRows && rowHdrs[row.rowid].used);
	return rowHdrs[row.rowid].rh;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline ColHdr& SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::colHdr(unsigned int colNo)
{
	assert(colNo < numCols && colHdrs[colHdrs[colNo].col].used);
	return colHdr(getCol(colNo));
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
inline ColHdr& SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::colHdr(Col col)
{
	assert(col.colid < maxNumCols && colHdrs[col.colid].colno < numCols && colHdrs[col.colid].used);
	return colHdrs[col.colid].ch;
}

/**
 * Prints the matrix in a (semi-)human-readable format to a given ostream.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::printMatrix(std::ostream& os) const
{
	os << "Matrix with " << numRows << " rows and " << numCols << " columns:" << std::endl;
	for(unsigned int i = 0; i < numRows; i++)
	{
		for(unsigned int j = 0; j < numCols; j++)
			os << std::setw(3) << (operator()(i, j) == PLUS ? "+1" : (operator()(i, j) == MINUS ? "-1" : "0"));
		os << std::endl;
	}
	os << std::endl;
}

/**
 * Checks whether two specified columns of the matrix
 * are equal (this doesn't compare the headers!).
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
bool SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::columnsEqual(Col col1, Col col2)
{
	assert(col1.colid < maxNumCols && col2.colid < maxNumCols && colHdrs[col1.colid].colno < numCols && colHdrs[col2.colid].colno < numCols && colHdrs[col1.colid].used && colHdrs[col2.colid].used);
	for(unsigned int i = 0; i < numRows; i++)
		if(operator()(getRow(i), col1) != operator()(getRow(i), col2))
			return false;
	return true;
}

/**
 * Writes zeros to a given column. If the column
 * with the specified index does not yet exist,
 * it is created. Returns the changed
 * or newly created column.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertZeroCol(unsigned int colNo, ColHdr colHdr)
{
	Col col = insertCol(colNo, colHdr);
	memset(matrix + col.colid * (maxNumRows / INTERNAL_SIZE), 0, (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
	return col;
}

/**
 * Writes the i-th unit vector to a given column.
 * If the column with the specified index does not
 * yet exist, it is created. Returns
 * the changed or newly created column.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertUnitCol(unsigned int colNo, unsigned int i, ColHdr colHdr)
{
	assert(i < numRows);
	Col col = insertCol(colNo, colHdr);
	memset(matrix + col.colid * (maxNumRows / INTERNAL_SIZE), 0, (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
	Row row; row.rowid = rowHdrs[i].row;
	setEntry(row, col, 1);
	return col;
}

/**
 * Inverts a specified column of the matrix and
 * writes the result to another (or the same) column. If the
 * latter doesn't exist, it is created. Returns
 * the changed or newly created column.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertInvCol(Col sourceCol, unsigned int targetColNo, ColHdr targetColHdr)
{
	assert(sourceCol.colid < maxNumCols && colHdrs[sourceCol.colid].colno < numCols && colHdrs[sourceCol.colid].used);
	Col col = insertCol(targetColNo, targetColHdr);
	for(unsigned int i = 0; i < (maxNumRows / INTERNAL_SIZE); i++)
	{
		unsigned int t = col.colid * (maxNumRows / INTERNAL_SIZE) + i;
		matrix[t] = ~matrix[sourceCol.colid * (maxNumRows / INTERNAL_SIZE) + i];
		INTERNAL_TYPE x = ((matrix[t] & mask1) << 1) & (matrix[t] & mask2);
		matrix[t] = matrix[t] ^ x ^ (x >> 1);
	}
	return col;
}

/**
 * Adds the values of two given columns and writes
 * it to another column. If the latter doesn't exist,
 * it is created and returned. Assumes that
 * the two columns can be added without leaving {-1,0,+1}
 * otherwise the behavior of this method is undefined.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertColSum(Col sourceCol1, Col sourceCol2, unsigned int targetColNo, ColHdr targetColHdr)
{
	assert(sourceCol1.colid < maxNumCols && sourceCol2.colid < maxNumCols && colHdrs[sourceCol1.colid].colno < numCols && colHdrs[sourceCol2.colid].colno < numCols && colHdrs[sourceCol1.colid].used && colHdrs[sourceCol2.colid].used);
	Col col = insertCol(targetColNo, targetColHdr);
	for(unsigned int i = 0; i < (maxNumRows / INTERNAL_SIZE); i++)
	{
		unsigned int t = col.colid * (maxNumRows / INTERNAL_SIZE) + i;
		matrix[t] = matrix[sourceCol1.colid * (maxNumRows / INTERNAL_SIZE) + i] | matrix[sourceCol2.colid * (maxNumRows / INTERNAL_SIZE) + i];
		INTERNAL_TYPE x = ((matrix[t] & mask1) << 1) & (matrix[t] & mask2);
		matrix[t] = matrix[t] ^ x ^ (x >> 1);
	}
	return col;
}

/**
 * Copies the values of a specified column (EXcluding
 * the header) to another (or the same) column. If the
 * latter doesn't exist, it is created. Returns the
 * index of the resulting column.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertColCopy(Col sourceCol, unsigned int targetColNo)
{
	assert(sourceCol.colid < maxNumCols && colHdrs[sourceCol.colid].colno < numCols && colHdrs[sourceCol.colid].used);
	Col col = insertCol(targetColNo, ColHdr());
	memcpy(matrix + col.colid * (maxNumRows / INTERNAL_SIZE), matrix + sourceCol.colid * (maxNumRows / INTERNAL_SIZE), (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
	return col;
}


template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertIntersectionCol(Col sourceCol1, Col sourceCol2, unsigned int targetColNo, ColHdr targetColHdr )
{
	assert(sourceCol1.colid < maxNumCols && sourceCol2.colid < maxNumCols && colHdrs[sourceCol1.colid].colno < numCols && colHdrs[sourceCol2.colid].colno < numCols && colHdrs[sourceCol1.colid].used && colHdrs[sourceCol2.colid].used);
	Col col = insertCol(targetColNo, targetColHdr);
	for(unsigned int i = 0; i < (maxNumRows / INTERNAL_SIZE); i++)
	{
		unsigned int t = col.colid * (maxNumRows / INTERNAL_SIZE) + i;
		matrix[t] = matrix[sourceCol1.colid * (maxNumRows / INTERNAL_SIZE) + i] & matrix[sourceCol2.colid * (maxNumRows / INTERNAL_SIZE) + i];
	}
	return col;
}


/**
 * Adds the entries of one column to another. Assumes that
 * the addition can be done without getting values outside
 * {-1,0,+1}.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::addToCol(Col sourceCol, Col targetCol)
{
	assert(sourceCol.colid < maxNumCols && targetCol.colid < maxNumCols && colHdrs[sourceCol.colid].colno < numCols && colHdrs[targetCol.colid].colno < numCols && colHdrs[sourceCol.colid].used && colHdrs[targetCol.colid].used);
	for(unsigned int i = 0; i < (maxNumRows / INTERNAL_SIZE); i++)
	{
		unsigned int t = targetCol.colid * (maxNumRows / INTERNAL_SIZE) + i;
		matrix[t] |= matrix[sourceCol.colid * (maxNumRows / INTERNAL_SIZE) + i];
		INTERNAL_TYPE x = ((matrix[t] & mask1) << 1) & (matrix[t] & mask2);
		matrix[t] = matrix[t] ^ x ^ (x >> 1);
	}
}

/**
 * Deletes a column from the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::deleteCol(Col col)
{
	assert(col.colid < maxNumCols && colHdrs[col.colid].colno < numCols && colHdrs[col.colid].used);
	// Decrease all column numbers greater that that of col by one.
	for(unsigned int i = colHdrs[col.colid].colno + 1; i < numCols; i++)
	{
		colHdrs[colHdrs[i].col].colno--;
		colHdrs[i - 1].col = colHdrs[i].col;
	}
	// Mark col as unused.
	colHdrs[col.colid].nextUnused = firstUnusedCol;
	firstUnusedCol = col.colid;
	colHdrs[col.colid].used = false;
	numCols--;
}

/**
 * Move a column including its header to a new
 * position in the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::moveCol(Col col, unsigned int newColNo)
{
	assert(col.colid < maxNumCols && colHdrs[col.colid].colno < numCols && newColNo < numCols);
	if(colHdrs[col.colid].colno > newColNo)
	{
		for(unsigned int i = colHdrs[col.colid].colno; i > newColNo; i--)
		{
			colHdrs[colHdrs[i - 1].col].colno++;
			colHdrs[i].col = colHdrs[i - 1].col;
		}
		colHdrs[col.colid].colno = newColNo;
		colHdrs[newColNo].col = col.colid;
	}
	else if(colHdrs[col.colid].colno < newColNo)
	{
		for(unsigned int i = colHdrs[col.colid].colno + 1; i <= newColNo; i++)
		{
			colHdrs[colHdrs[i].col].colno--;
			colHdrs[i - 1].col = colHdrs[i].col;
		}
		colHdrs[col.colid].colno = newColNo;
		colHdrs[newColNo].col = col.colid;
	}
	// else: nothing to be done
}

/**
 * Fills a specified row with zeros after creating
 * it if necessary.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Row SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertZeroRow(unsigned int rowNo, RowHdr rowHdr)
{
	Row row = insertRow(rowNo, rowHdr);
	for(unsigned int i = 0; i < numCols; i++)
		setEntry(row, getCol(i), 0);
	return row;
}

/**
 * Inserts a copy of a given row (including its
 * header) at another position in the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Row SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertRowCopy(Row sourceRow, unsigned int targetRowNo)
{
	assert(sourceRow.rowid < maxNumRows && rowHdrs[sourceRow.rowid].rowno < numRows && rowHdrs[sourceRow.rowid].used);
	Row row = insertRow(targetRowNo, rowHdrs[sourceRow.rowid].rh);
	for(unsigned int i = 0; i < numCols; i++)
		setEntry(row, getCol(i), operator()(sourceRow, getCol(i)));
	return row;
}

/**
 * Deletes a row from the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::deleteRow(Row row)
{
	assert(row.rowid < maxNumRows && rowHdrs[row.rowid].rowno < numRows && rowHdrs[row.rowid].used);
	// Decrease all row numbers greater that that of row by one.
	for(unsigned int i = rowHdrs[row.rowid].rowno + 1; i < numRows; i++)
	{
		rowHdrs[rowHdrs[i].row].rowno--;
		rowHdrs[i - 1].row = rowHdrs[i].row;
	}
	// Mark row as unused.
	rowHdrs[row.rowid].nextUnused = firstUnusedRow;
	firstUnusedRow = row.rowid;
	rowHdrs[row.rowid].used = false;
	numRows--;
}

/**
 * Move a row including its header to a new
 * position in the matrix.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::moveRow(Row row, unsigned int newRowNo)
{
	assert(row.rowid < maxNumRows && rowHdrs[row.rowid].rowno < numRows && newRowNo < numRows);
	if(rowHdrs[row.rowid].rowno > newRowNo)
	{
		for(unsigned int i = rowHdrs[row.rowid].rowno; i > newRowNo; i--)
		{
			rowHdrs[rowHdrs[i - 1].row].rowno++;
			rowHdrs[i].row = rowHdrs[i - 1].row;
		}
		rowHdrs[row.rowid].rowno = newRowNo;
		rowHdrs[newRowNo].row = row.rowid;
	}
	else if(rowHdrs[row.rowid].rowno < newRowNo)
	{
		for(unsigned int i = rowHdrs[row.rowid].rowno + 1; i <= newRowNo; i++)
		{
			rowHdrs[rowHdrs[i].row].rowno--;
			rowHdrs[i - 1].row = rowHdrs[i].row;
		}
		rowHdrs[row.rowid].rowno = newRowNo;
		rowHdrs[newRowNo].row = row.rowid;
	}
	// else: nothing to be done
}

/**
 * Produces a copy of the matrix using only the given columns.
 * The row and column headers are NOT copied.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::clone(const SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>& mat, std::list<Col>& colList)
{
	// Delete this matrix
	free();
	// Determine size
	numRows = mat.numRows;
	numCols = colList.size();
 	maxNumRows = mat.maxNumRows;
 	maxNumCols = numCols;
 	// Copy row headers
 	rowHdrs.resize(maxNumRows);

 	firstUnusedRow = mat.firstUnusedRow;
 	for(unsigned int i = 0; i < maxNumRows/*numRows*/; i++)
 	 	//	rowHdrs[i] = mat.rowHdrs[i];
 	{
 		rowHdrs[i].row = mat.rowHdrs[i].row;
 		rowHdrs[i].rowno = mat.rowHdrs[i].rowno;
 		rowHdrs[i].used = mat.rowHdrs[i].used;
 		rowHdrs[i].nextUnused = mat.rowHdrs[i].nextUnused;
 	}

 	// Column headers are copied later, for now just make room for them
 	colHdrs.resize(maxNumCols);
 	if(numCols < maxNumCols)
 		firstUnusedCol = numCols;
 	else
 		firstUnusedCol = -1;
 	for(unsigned int i = numCols; i < maxNumCols; i++)
 	{
 		colHdrs[i].used = false;
		colHdrs[i].nextUnused = i + 1;
	}
	colHdrs[maxNumCols-1].nextUnused = -1;
	// Copy matrix
 	if(maxNumRows * maxNumCols > 0)
 	{
 		matrix = new INTERNAL_TYPE[maxNumCols * (maxNumRows / INTERNAL_SIZE)];
 		assert(matrix != NULL);
 		memset(matrix, 0, maxNumCols * (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
 		unsigned int i = 0;
 		for(std::list<Col>::iterator iter = colList.begin(); iter != colList.end(); iter++)
 		{
 			// Copy column
 			memcpy(matrix + i * (maxNumRows / INTERNAL_SIZE), mat.matrix + iter->colid * (mat.maxNumRows / INTERNAL_SIZE), (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
 			// Set column header
 			colHdrs[i].used = true;
 			//colHdrs[i].ch = mat.colHdr(*iter);
 			colHdrs[i].colno = i;
 			colHdrs[i].col = i;
 			iter->colid = i;
 			i++;
 		}
 	}
}

/**
 * Copies a given matrix into this one as a block starting
 * at the given offset.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertMatrix(const SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>& mat, unsigned int rowOffset, unsigned int colOffset)
{
	//FIXME: This needs to be sped up!
	throw " should not use this function insertMatrix!";
	if(rowOffset + mat.numRows < maxNumRows || colOffset + mat.numCols < maxNumCols)
		extend(rowOffset + mat.numRows, colOffset + mat.numCols);

	for(unsigned int c = 0; c < mat.numCols; c++)
		insertZeroCol(colOffset + c, mat.colHdrs[mat.colHdrs[c].col].ch);
	for(unsigned int r = 0; r < mat.numRows; r++)
		insertZeroRow(rowOffset + r, mat.rowHdrs[mat.rowHdrs[r].row].rh);
	for(unsigned int c = 0; c < mat.numCols; c++)
		for(unsigned int r = 0; r < mat.numRows; r++)
			setEntry(rowOffset + r, colOffset + c, mat.operator()(r, c));
}


/**
 * Extends the matrix to have at least the given number of
 * rows and columns.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
void SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::extend(unsigned int rows, unsigned int cols)
{
	unsigned int	newNumRows, newNumCols;
	if(rows * EXTEND_CONDITION > maxNumRows)
	{
		newNumRows = std::max(rows, (unsigned int)ceil((maxNumRows * EXTEND_FACTOR)));
		newNumRows = ((newNumRows + INTERNAL_SIZE - 1) / INTERNAL_SIZE) * INTERNAL_SIZE;
	}
	else
		newNumRows = maxNumRows;
	if(cols * EXTEND_CONDITION > maxNumCols)
		newNumCols = std::max(cols, (unsigned int)ceil((maxNumCols * EXTEND_FACTOR)));
	else
		newNumCols = maxNumCols;
	// Reallocate
	if(newNumCols > maxNumCols || newNumRows > maxNumRows)
	{
		INTERNAL_TYPE* newMatrix = NULL;
		if((unsigned long long)newNumCols * (unsigned long long)newNumRows > 0)
		{
			newMatrix = new INTERNAL_TYPE[newNumCols * (newNumRows / INTERNAL_SIZE)];
			memset(newMatrix, 0, newNumCols * (newNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
			for(unsigned int j = 0; j < numCols; j++)
				//for(unsigned int i = 0; i < numRows; i++)
				//	newMatrix[i + j * newNumRows] = operator()(i, j);
				memcpy(newMatrix + j * (newNumRows / INTERNAL_SIZE), matrix + j * (maxNumRows / INTERNAL_SIZE), (maxNumRows / INTERNAL_SIZE) * sizeof(INTERNAL_TYPE));
		}
		if(matrix != NULL)
			delete [] matrix;
		matrix = newMatrix;
	}
	if(newNumRows > maxNumRows)
	{
		rowHdrs.resize(newNumRows);

		for(unsigned int i = maxNumRows; i < newNumRows; i++)
		{
			rowHdrs[i].used = false;
			rowHdrs[i].nextUnused = i + 1;
		}
		rowHdrs[newNumRows-1].nextUnused = firstUnusedRow;
		firstUnusedRow = maxNumRows;
		maxNumRows = newNumRows;
	}
	if(newNumCols > maxNumCols)
	{
		colHdrs.resize(newNumCols);
		for(unsigned int i = maxNumCols; i < newNumCols; i++)
		{
			colHdrs[i].used = false;
			colHdrs[i].nextUnused = i + 1;
		}
		colHdrs[newNumCols-1].nextUnused = firstUnusedCol;
		firstUnusedCol = maxNumCols;

		maxNumCols = newNumCols;
	}
}

/**
 * Internal methods whose job it is to find or create
 * unused space for a new row or column and create it.
 * It also adjusts the header to make it appear as if
 * the new row or column was inserted at the given
 * position.
 */
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Row SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertRow(unsigned int rowNo, RowHdr rowHdr)
{
	assert(rowNo <= numRows);
	//checkExtend(numRows + 1, numCols);
	if(firstUnusedRow == -1)
	{
		assert(numRows == maxNumRows);
		extend(numRows + 1, numCols);
	}
	assert(firstUnusedRow != -1);

	int n = firstUnusedRow;
	firstUnusedRow = rowHdrs[n].nextUnused;

	assert(!rowHdrs[n].used);

	// Increment all row numbers >= rowNo by one
	for(unsigned int j = numRows; j > rowNo; j--)
	{
		rowHdrs[j].row = rowHdrs[j - 1].row;
		rowHdrs[rowHdrs[j].row].rowno++;
	}
	// Insert the new row
	rowHdrs[n].used = true;
	rowHdrs[n].rowno = rowNo;
	rowHdrs[n].rh = rowHdr;
	rowHdrs[rowNo].row = n;
	numRows++;
	Row row; row.rowid = n;
	return row;
}
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
Col SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::insertCol(unsigned int colNo, ColHdr colHdr)
{
	assert(colNo <= numCols);

	if(firstUnusedCol == -1)
	{
		assert(numCols == maxNumCols);
		extend(numRows , numCols+ 1);
	}
	assert(firstUnusedCol != -1);

	int n = firstUnusedCol;
	firstUnusedCol = colHdrs[n].nextUnused;

	assert(!colHdrs[n].used);
	// Increment all row numbers >= rowNo by one
	for(unsigned int j = numCols; j > colNo; j--)
	{
		colHdrs[j].col = colHdrs[j - 1].col;
		colHdrs[colHdrs[j].col].colno++;
	}

	// Insert the new column
	colHdrs[n].used = true;
	colHdrs[n].colno = colNo;
	colHdrs[n].ch = colHdr;
	colHdrs[colNo].col = n;
	numCols++;
	Col col; col.colid = n;
	return col;
}

#ifdef DEBUG
template<class RowHdr, class ColHdr, typename INTERNAL_TYPE>
bool SignMatrix<RowHdr, ColHdr, INTERNAL_TYPE>::checkInternalIntegrity()
{
	for(unsigned int i = 0; i < numRows; i++)
	{
		for(unsigned int j = 0; j < numRows; j++)
		{
			if(operator()(i, j) != ZERO && operator()(i, j) != PLUS && operator()(i, j) != MINUS)
			{
				std::cout << "Matrixfehler in Zeile " << i << ", Spalte " << j << std::endl;
				std::cout.flush();
				assert(false);
			}
		}
	}
}
#endif

}

#endif /* SIGNMATRIX_H_ */
