import Select from 'react-select'

export const RenameItem = props => {

  // Get from domain options based on domains in the rule
  const options = props.rule.domains.map(d => ({ label: d, value: d }))

  const handleChangeFrom = event => {
    props.handleChange({
      target: { name: "from", value: event.value }
    })
  }
  const handleChangeAfter = event => {
    props.handleChange({
      target: { name: "after", value: event ? event.map(e => e.value) : [] }
    })
  }

  return (
    <li>
      <button
        type="button"
        onClick={props.handleRemove}
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
          value={options.find(o => o.label === props.data.from)}
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
          value={props.data.after.map(a => options.find(o => o.label === a))}
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
          value={props.data.to}
          onChange={props.handleChange}
        />
      </div>
    </li>
  )
}
