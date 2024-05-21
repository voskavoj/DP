"""
    Text visualization utilities
"""

def tableify(data, col_dec=None, col_unit=None, col_head=None, row_head=None):
    try:
        rows = list()

        if col_head:
            column_header = ""
            for ch in col_head:
                column_header += f"{ch} & "
            column_header = column_header[:-3] + "\\\\ \\hline"
            rows.append(column_header)

        for j, row_data in enumerate(data):
            row = ""

            if row_head:
                row += f"{row_head[j]} & "

            for i, col in enumerate(row_data):
                if col_dec:
                    try:
                        if col_dec[i] == 0:
                            col = round(col)
                        else:
                            col = round(col, col_dec[i])
                    except Exception:
                        pass

                if col_unit:
                    cu = " " + col_unit[i]
                else:
                    cu = ""

                row += f"{col}{cu} & "
            row = row[:-3] + "\\\\"
            rows.append(row)

        print()
        for r in rows:
            print(r.replace("%", "\%").replace("_", "\_"))
        print()

    except Exception as e:
        print(f"Error in tableify: {e}")
