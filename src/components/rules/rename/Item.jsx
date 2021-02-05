import React, { useMemo } from 'react'
import Select from 'react-select'


const RenameItem = ({
  from,
  before,
  after,
  to,
  options,
  handleChange,
  handleRemove,
}) => {

  console.log(from, after, to, options.find(o => o.label === from))

  const handleChangeFrom = event => {
    handleChange({
      target: { name: "from", value: event.innerValue }
    })
  }
  const handleChangeAfter = event => {
    handleChange({
      target: { name: "after", value: event ? event.map(e => e.innerValue) : [] }
    })
  }
  const handleChangeBefore = event => {
    handleChange({
      target: { name: "before", value: event ? event.map(e => e.innerValue) : [] }
    })
  }

  const nameValue = useMemo(
    () => options ? options.find(o => o.label === from) : null,
    [from, options]
  )
  const afterValue = useMemo(
    () => options ? after.map(a => options.find(o => o.label === a)) : [],
    [after, options]
  )
  const beforeValue = useMemo(
    () => options ? before.map(b => options.find(o => o.label === b)) : [],
    [before, options]
  )

  return (
    <li>
      <button
        type="button"
        onClick={handleRemove}
      >
        Delete
      </button>

      {/* Change this domain name */}
      <div className="rule-field">
        <label htmlFor="renameName">From:</label>
        <Select
          id="renameName"
          className="select"
          onChange={handleChangeFrom}
          options={options}
          value={nameValue}
        />
      </div>

      {/* Change domain name when occuring before these domains */}
      <div className="rule-field">
        <label htmlFor="renameBefore">Before domains:</label>
        <Select
          id="renameBefore"
          className="select"
          options={options}
          onChange={handleChangeBefore}
          value={beforeValue}
          isMulti
        />
      </div>

      {/* Change domain name when occuring after these domains */}
      <div className="rule-field">
        <label htmlFor="renameAfter">After domains:</label>
        <Select
          id="renameAfter"
          className="select"
          options={options}
          onChange={handleChangeAfter}
          value={afterValue}
          isMulti
        />
      </div>

      {/* Change domain name to this value */}
      <div className="rule-field">
        <label htmlFor="filterTo">To:</label>
        <input
          id="filterTo"
          type="text"
          name="to"
          value={to}
          onChange={handleChange}
        />
      </div>
    </li>
  )
}

export default React.memo(RenameItem)
