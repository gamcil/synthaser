import React from 'react'

import RenameItem from './Item'

const RenameList = ({
  renames,
  options,
  handleAdd,
  handleChange,
  handleRemove,
}) => (
  <div>
    <div>
      <button type="button" onClick={handleAdd}>Add</button>
    </div>
    <ul>
      {renames.map(rename => (
        <RenameItem
          key={rename.uuid}
          from={rename.from}
          before={rename.before}
          after={rename.after}
          to={rename.to}
          options={options}
          handleChange={handleChange(rename.uuid)}
          handleRemove={handleRemove(rename.uuid)}
        />
      ))}
    </ul>
  </div>
)

export default React.memo(RenameList)
